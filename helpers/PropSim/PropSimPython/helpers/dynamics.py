# Dynamics modelling methods for PropSimPython
# Includes 
# Project Caelus, Aphlex 1C Engine
# Jason Chen, 20 February, 2021


import numpy as np
from classes import Struct
from nozzle import nozzle_calc
from n2o import n2o_properties
from main_functs import find_G


def n2o_tank_mdot(inputs: Struct, state: LiquidStateVector, time: float) -> tuple:
    """
    Calculates flow rate out of N2O Tank. Interpolates between liquid and gas flow 
    based on amount of liquid oxidizer remaining. If no liquid oxidizer remains, it
    uses gas flow. If plenty of liquid oxidizer (i.e. mass greater than set tolerance),
    liquid flow is used. Linear interpolation is used (in order to avoid hysteresis)
    in between with a small, steady liquid source in tank. 

    Inputs:
        - inputs: object representing motor inputs
        - state: object representing state of system
        - time: time
    Outputs:
        - m_dot_lox: mass flow rate out of tank of liquid oxidizer
        - m_dot_gox: mass flow rate out of tank of gaseous oxidizer
        - m_dot_oxtank_press: mass flow rate out of tank of pressurant gas
        - T_dot_drain: oxidizer tank temperature change from draining
        - p_crit: critical pressure below which draining flow is choked
        - m_dot_ox_crit: critical mass flow rate below which draining flow is choked
    """
    # Constants
    dm_lox_tol = 1e-2 # The tolerance as mentioned above for "plenty" of liquid [kg]
    M_n2o = 0.044013 # Molecular mass of nitrous oxide [kg/mol]
    a_n2o = 0.38828/M_n2o^2 # van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]

    inj_A_eff = inputs.ox.injector_area*inputs.throttle(time) # Effective area of the injector

    # If tank pressure is greater than combustion chamber pressure
    if state.p_oxmanifold > state.p_cc and inputs.ox.Cd_injector*inj_A_eff > 0:
        # Find flow rates for gas, flow rates for liquid, interpolate within tolerance for smooth transition
        m_dot_lox          = np.zeros((2,1))
        m_dot_gox          = np.zeros((2,1))
        m_dot_oxtank_press = np.zeros((2,1))
        p_crit             = np.zeros((2,1))
        m_dot_ox_crit      = np.zeros((2,1))
        Q                  = np.zeros((2,1)) # Volumetric flow rate

        # Liquid flow
        m_dot_ox, m_dot_ox_crit[0], p_crit[0] = liq_n2o_mdot(inputs, inj_A_eff, 
            state.p_oxmanifold, state.T_oxtank, state.p_cc)

        m_dot_lox[0] = m_dot_ox
        m_dot_gox[0] = 0
        m_dot_oxtank_press[1] = 0
        Q[1] = m_dot_ox/state.n2o_props.rho_l
        
        # Gas Flow
        d_inj = np.sqrt(4/np.pi*inj_A_eff)
        _,  _,  _,  _, m_dot_ox, _ = nozzle_calc( d_inj, d_inj, state.T_oxtank, \
            state.p_oxmanifold, state.gamma_ox_ullage, M_n2o, state.p_cc)
        
        p_crit[1] = 0
        m_dot_ox_crit[1] = 0
        m_dot_lox[1] = 0
        m_dot_gox[1] = m_dot_ox*\
            (state.m_gox)/(state.m_gox + state.m_oxtank_press)
        m_dot_oxtank_press[1] = m_dot_ox*\
            (state.m_oxtank_press)/(state.m_gox + state.m_oxtank_press)
        Q[1] = m_dot_ox*state.V_ox_ullage/(state.m_gox + state.m_oxtank_press)
        
        # Total Flow Rate
        if state.m_lox > dm_lox_tol:
            frac_lox = 1
        else:
            frac_lox = max(0,state.m_lox/dm_lox_tol)
        
        m_dot_lox = frac_lox*m_dot_lox(1) + (1-frac_lox)*m_dot_lox[1]
        m_dot_gox = frac_lox*m_dot_gox(1) + (1-frac_lox)*m_dot_gox[1]
        m_dot_oxtank_press = frac_lox*m_dot_oxtank_press(1) + (1-frac_lox)*m_dot_oxtank_press[1]
        p_crit = frac_lox*p_crit(1) + (1-frac_lox)*p_crit[1]
        m_dot_ox_crit = frac_lox*m_dot_ox_crit(1) + (1-frac_lox)*m_dot_ox_crit[1]
        Q = frac_lox*Q[0] + (1-frac_lox)*Q[1]
    else: # bruh ur bad if ur chamber pressure is greater than tank press
        m_dot_lox = 0
        m_dot_gox = 0
        m_dot_oxtank_press = 0
        p_crit = state.p_oxtank
        m_dot_ox_crit = 0
        Q = 0


def liq_n2o_mdot(inputs: Struct, A_inj: float, P_up: float, T_up: float, P_down: float) -> tuple:
    """ 
    Calculation of the mass flow rate of liquid nitrous oxide.
    Inputs:
        - inputs: object representing motor inputs
        - A_inj: effective injector area
        - P_up: pressure upstream
        - T_up: temperature upstream
        - P_down: pressure downstream
    Outputs:
        - m_dot_ox: mass flow rate of the oxidizer (N2O)
        - m_dot_crit: critical mass flow rate of the oxidizer (N2O)
        - p_crit: critical pressure (between gas and liquid form)
    """
    # Set pressure to be a minimum of vapor pressure to address modeling limitation.
    temp_n2o_props = n2o_properties(T_up)
    if P_up < temp_n2o_props.p_vap:
        P_up = temp_n2o_props.p_vap
    
    # Use 2-phase model for low supercharge pressure ratio
    if P_up/temp_n2o_props.Pvap < 3:
        G, G_crit, p_down_crit = find_G(T_up, P_up, P_down)
    else:
        G = np.sqrt(2*(P_up-P_down)*temp_n2o_props.rho_l)
        G_crit = 0
        p_down_crit = 0

    # Find the effective CdA throughout the system
    # Cv_mpv = 0.7 % MPV value of the main propellant valve
    # CdA_mpv = Cv_mpv*7.598e-7/sqrt(2/1000)  # 1000 kg/m^3, density of water
    # injector_area = inputs.ox.Cd_injector*A_inj
    # CdA_eff = 1/sqrt((1/(CdA_mpv^2))+(1/(injector_area^2)))

    # m_dot_ox = CdA_eff*G # inputs.ox.Cd_injector*A_inj*G
    # m_dot_crit = CdA_eff*G_crit
    m_dot_ox = inputs.ox.Cd_injector*A_inj*G
    m_dot_crit = inputs.ox.Cd_injector*A_inj*G_crit
    p_crit = p_down_crit
