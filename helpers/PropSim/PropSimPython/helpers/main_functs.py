# Main helper methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# Jason Chen, 10 February, 2021


import numpy as np
from scipy.integrate import solve_ivp
from typing import Tuple
from classes import Struct
from state_flow import init_liquid_state, LiquidStateVector, lsv_from_column_vec
from dynamics import n2o_tank_mdot
from nozzle import nozzle_calc
from design_liquid import design_liquid
from masscalc import calc_mass


def integration(inputs: Struct) -> Tuple[float, Struct]:
    """
    Integrate necesary differential equations for rocket engine modeling using Euler's method. 
    TODO: Include real combustion properties and supercharging.
    INPUTS:
        - inputs: structure of motor characteristics (all units SI: m, s, 
        kg, K, mol)
            - CombustionData: string filename in "Combustion Data/" folder 
            from which to source combustion data
            - Injexit_area: total orifice area of injector
            - Cdischarge: discharge coefficient of injector orifices
            - MOVTime: time for injector flow to ramp up to 100
            - rocket_dry_mass: dry mass of rocket
            - tankvol: volume of oxidizer tank
            - l_vol: initial volume of liquid nitrous oxide in oxidizer 
            tank
            - tank_id: inner diameter of tank
            - h_offset: height difference between bottom of tank and 
            injector
            - flowline_id: inner diameter of flow line
            - Ttank: tank initial temperature
            - Fueldensity: density of fuel
            - grainlength: length of fuel grain
            - graindiameter: outer diameter of grain
            - portradius0: grain intial port radius
            - chamberlength: combustion chamber total length
            - n: grain ballistic coefficient (r_dot = a*G^n)
                - a: grain ballistic coefficient (r_dot = a*G^n)
            - nozzle_efficiency: nozzle exhaust energy efficiency 
            (proportion)
            - nozzle_correction_factor: proportion of ideal thrust 
            (factor for divergence losses, etc.) actually achieved
            - c_star_efficiency: combustion efficiency / proportion of c*
            actually achieved
            - Tdiameter: nozzle throat diameter
            - E: expansion ratio
            - Tamb: ambient temperature
            - Pamb: ambient pressure
            - SPress: supercharging regulator pressure
            - M_sc: molecular mass
            - SVol: volume of external pressurization tank (0 for 
            supercharging / no external tank)
            - c_v_S: specific heat at constant volume of pressurant gas
                - c_p_S: specific heat at constant volume of pressurant gas
            - P_S_i: initial pressurant gas storage pressure (must be 
            present, but not used for supercharging / no external tank)
            - S_CdA: effective flow area of pressurant gas (must be 
            but not used for supercharging / no external tank)
            - Charged: 1 for pressurant gas present, 0 for no pressurant 
            gas present
        - mode: structure of options defining mode of motor operation
            - combustion_on: 1 for hot-fire, 0 for cold-flow
            - flight_on: 1 for flight conditions, 0 for ground conditions
        - tspan: vector of time values over which to record outputs
    OUTPUS:
        - tspan: output time vector
        - F_thrust: thrust over tspan
        - p_cc: combustion chamber pressure over tspan
        - p_oxtank: tank pressure over tspan
        - p_oxmanifold: oxidizer manifold pressure over tspan
        - T_tank: tank temperature over tspan
        - T_cc: combustion chamber as temperature over tspan
        - area_core: fuel grain core area over tspan
        - OF: oxidizer/fuel ratio over tspan
        - m_dot_ox: oxidizer mass flow rate over tspan
        - m_dot_ox_crit: critical two-phase oxidizer mass flow rate over tspan
        - p_crit: two-phase critical downstream pressure over tspan
        - M_e: exit Mach number over tspan
        - p_exit: nozzle exit pressure over tspan
        - p_shock: critical back pressure for normal shock formation over
        tspan
    """
    # Recording the output data for this timestep
    record_config = { 
        "F_thrust"        : None,
        "p_cc"            : None,
        "p_oxtank"        : None,
        "p_oxpresstank"   : None,
        "p_fueltank"      : None,
        "p_fuelpresstank" : None,
        "p_oxmanifold"    : None,
        "T_oxtank"        : None,
        "T_cc"            : None,
        "area_core"       : None,
        "of_ratio_i"      : None,
        "gamma_ex"        : None,
        "m_dot_ox"        : None,
        "m_dot_fuel"      : None,
        "p_crit"          : None,
        "m_dot_ox_crit"   : None,
        "M_e"             : None,
        "p_exit"          : None,
        "p_shock"         : None
    } 
    record = Struct(record_config) # Convert into Struct class
    state_0, x0 = init_liquid_state(inputs)
    # Configure the integrator 
    solve_ivp(liquid_model, method='BDF') # NOTE: Incomplete method call
    

def liquid_model(time: float, x: np.ndarray, inputs: Struct):
    """ Model engine physics for a liquid motor, provide state vector derivative. """
    
    # Constants
    R_u = inputs.constants.R_u # Universal gas constant [J/mol*K]
    M_n2o = 0.044013 # Molecular mass of nitrous oxide [kg/mol]
    R_n2o = R_u/M_n2o #  Specific gas constant of nitrous oxide [J/kg*K]
    a_n2o = 0.38828/M_n2o^2 # van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]
    b_n2o = 44.15/M_n2o*10^-6 # van der Waal's constant a for N2O [m^3/kg]

    state = lsv_from_column_vec(x, inputs)
    state_dot = LiquidStateVector()

    # Calculate injector mass flow rate
    # TODO: This is where Python conversion stopped (02/25/2021); see n2o.py for n2o_tank_mdot() work
    output = n2o_tank_mdot(inputs, state, time)
    m_dot_lox, m_dot_gox, m_dot_oxtank_press, T_dot_drain_ox, p_crit, m_dot_ox_crit = output


def find_G():
    raise NotImplementedError


# TODO: The following functions (find_mass_gradient and create_struct) are UNTESTED
def find_mass_gradient(constraints, values: np.ndarray, delta, results_prev):
    goal = create_struct(constraints, values)[0]
    design = create_struct(constraints, values)[1]

    base_inputs = design_liquid(results_prev, goal, design, False)
    central_value = calc_mass(base_inputs)

    grad = np.zeros((1, len(values)))
    
    for j in range(len(values)):
        step = np.zeros(np.ndarray.ndim(values))
        step[j] = values[j]*delta

        goal = create_struct(constraints, values + step)[0]
        design = create_struct(constraints, values + step)[1]
        grad[j] = (calc_mass(design_liquid(base_inputs, goal, design, False)) - central_value)/delta
    grad = np.divide(grad,np.linalg.norm(grad))
    return grad


def create_struct(constraints, values: np.ndarray):
    goal['total_impulse'] = constraints[1]
    goal['max_thrust'] = constraints[2]
    goal['OF'] = constraints[3]
    goal['min_fuel_dp'] = constraints[4]
    goal['min_ox_to_dp'] = constraints[5]
    goal['ox_to_fuel_time'] = values[1]
    design['ox_ullage'] = constraints[6]
    design['p_tanks'] = values[2]
    design['exp_ratio'] = values[3]
    goal = Struct(goal)
    design = Struct(design)
    return goal, design
