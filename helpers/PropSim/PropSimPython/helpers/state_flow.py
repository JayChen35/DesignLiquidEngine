# State-flow helper methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# Jason Chen, 10 February, 2021


import numpy as np
from classes import Struct
from n2o import n2o_properties

class StateVector():
    """
    A class representing an engine state vector, allowing access to variables by name.
    """
    def __init__(self, inputs: Struct):
        super().__init__()
        self.inputs = inputs
        # Tank properties
        self.m_lox = None # Mass of liquid oxidizer
        self.m_gox = None # Mass of gaseous Oxidizer
        self.m_oxtank_press = None # Mass of pressurant in oxidizer tank
        self.m_oxpresstank = None # Mass of pressurant in oxidizer pressurant tank
        self.T_oxtank = None # Temperature oxidizer tank
        self.n2o_props = None # Nitrous oxide properties
        # Combustion chamber
        self.m_cc = None # Mass of gas (exhaust product in chamber)
        self.M_cc = None # Molecular mass of gas
        self.gamma_cc = None  # Ratio of specific heats of gas
        self.T_cc = None # Temperature of gas

    def V_lox(self):
        # Calculate volume of liquid oxidizer
        V_lox = self.m_lox / self.n2o_props.rho_l
        return V_lox

    def p_oxtank(self):
        # Calculate volume of liquid oxidizer
        p_oxtank = self.p_gox + self.p_oxtank_press
        return p_oxtank
        
    def oxtank_m_cv(self):
        # Thermal capacity of oxidizer tank (constant volume)
        m_cv = (self.m_lox*self.n2o_props.cv_l + self.m_gox*self.n2o_props.cv_g + \
            self.m_oxtank_press*self.inputs.ox_pressurant.gas_props.c_v)
        return m_cv
    
    def oxtank_m_cp(self):
        # Thermal capacity of oxidizer tank (constant pressure)
        m_cp = (self.m_lox*self.n2o_props.cp_l + self.m_gox*self.n2o_props.cp_g + \
            self.m_oxtank_press*self.inputs.ox_pressurant.gas_props.c_p)
        return m_cp
    
    def p_gox(self):
        # Calculate oxidizer tank partial pressure of oxidizer
        M_n2o = 0.044013 # Molecular mass of nitrous oxide [kg/mol]
        R_n2o = self.inputs.constants.R_u/M_n2o #  Specific gas constant of nitrous oxide [J/kg*K]
        a_n2o = 0.38828/M_n2o^2 # van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]
        b_n2o = 44.15/M_n2o*10^-6 # van der Waal's constant a for N2O [m^3/kg]
        p_gox = press_VDW(self.T_oxtank, self.m_gox/self.V_ox_ullage, R_n2o, a_n2o, b_n2o)
        return p_gox
    
    def p_oxmanifold(self):
        # Calculate manifold pressure
        accel_net = 0 # Net acceleration initially is zero, and then is modified later in integration.py
        # Calculate height difference between liquid level and manifold
        dh_lox = self.inputs.ox.h_offset_tank + self.V_lox/(np.pi*(self.inputs.ox.tank_id/2)**2)
        # Calculate acceleration of oxidizer
        if self.inputs.options.flight_on:
            accel_rel = accel_net + self.inputs.constants.g_0
        else:
            accel_rel = self.inputs.constants.g_0
        # Calculate manifold pressure
        p_oxmanifold = self.p_oxtank + accel_rel*self.n2o_props.rho_l*dh_lox
        return p_oxmanifold
    
    def p_oxtank_press(self):
        # Calculate oxidizer tank partial pressure of pressurant
        p_oxtank_press = self.m_oxtank_press*self.inputs.ox_pressurant.gas_props.R_spec* \
            self.T_oxtank/self.V_ox_ullage
        return p_oxtank_press
    
    def p_oxpresstank(self):
        # Calculate pressure in oxidizer pressurant tank
        m_press_initial = self.inputs.ox_pressurant.storage_init_press* \
            self.inputs.ox_pressurant.tank_V/( \
            self.inputs.ox_pressurant.gas_props.R_spec* \
            self.inputs.T_amb)
        p_oxpresstank = self.inputs.ox_pressurant.storage_init_press* \
            (self.m_oxpresstank/m_press_initial)^ \
            self.inputs.ox_pressurant.gas_props.gamma
        return p_oxpresstank
    
    def T_oxpresstank(self):
        # Calculate temperature in oxidizer pressurant tank
        T_oxpresstank = self.inputs.T_amb* \
            (self.p_oxpresstank/self.inputs.ox_pressurant.storage_init_press)^ \
            ((self.inputs.ox_pressurant.gas_props.gamma-1)/ \
            self.inputs.ox_pressurant.gas_props.gamma)
        return T_oxpresstank
    
    def V_ox_ullage(self):
        # Calculate ullage volume in oxidizer tank
        V_ox_ullage = self.inputs.ox.V_tank - self.V_lox
        return V_ox_ullage
    
    def gamma_ox_ullage(self):
        # Calculate ullage gas ratio of specific heats
        gamma_ox_ullage = (self.m_gox*self.n2o_props.cp_g +  \
            self.m_oxtank_press*self.inputs.ox_pressurant.gas_props.c_p)/ \
            (self.m_gox*self.n2o_props.cv_g +  \
            self.m_oxtank_press*self.inputs.ox_pressurant.gas_props.c_v)
        return gamma_ox_ullage
    
    def p_cc(self):
        # Calculate pressure in the combustion chamber
        R_u = 8.3144621 # Universal gas constant [J/mol*K]
        p_cc = (self.m_cc/self.V_cc)*(R_u/self.M_cc)*self.T_cc
        return p_cc

    def press_VDW(self, T: float, rho: float, R: float, a: float, b: float):
        # VanDerWaal's equation of state, solving for pressure
        p = R*T/(1/rho-b)-a*rho**2 # ./ is item-wise division simply / in NumPy.
        return p


class LiquidStateVector(StateVector):
    def __init__(self, inputs: dict):
        super().__init__()
        self.T_fueltank_press = None # Fuel tank pressurant temperature
        self.m_fueltank_press = None # Mass of pressurant in fuel tank
        self.m_fuelpresstank_press = None # Mass of pressurant in fuel pressurant tank
        self.m_fuel = None # Fuel mass, kg
        self.inputs = inputs
    
    def p_fuel_tank(self):
        # Calculate fuel tank pressure
        p_fueltank = self.m_fuel_tank_press*self.inputs.fuel_pressurant.gas_props.R_spec* \
            self.T_fueltank_press/self.V_fuel_ullage
        return p_fueltank
    
    def p_fuelpresstank(self):
        # Calculate pressure in fuel pressurant tank
        m_press_initial = self.inputs.fuel_pressurant.storage_init_press* \
            self.inputs.fuel_pressurant.tank_volume/ \
                (self.inputs.fuel_pressurant.gas_props.R_spec*self.inputs.T_amb)
        p_fuelpresstank = self.inputs.fuel_pressurant.storage_init_press* \
            (self.m_fuelpresstank/m_press_initial)^self.inputs.fuel_pressurant.gas_props.gamma
        return p_fuelpresstank
    
    def T_fuelpresstank(self):
        # Calculate temperature in oxidizer pressurant tank
        T_fuelpresstank = self.inputs.T_amb* \
            (self.p_fuelpresstank/self.inputs.fuel_pressurant.storage_init_press)** \
            ((self.inputs.fuel_pressurant.gas_props.gamma-1)/self.inputs.fuel_pressurant.gas_props.gamma)
        return T_fuelpresstank

    def p_fuelpresstank(self):
        p_fuelpresstank = 0
        return p_fuelpresstank

    def T_fuelpresstank(self):
        T_fuelpresstank = 298
        return T_fuelpresstank
    
    def V_fuel(self):
        # Calculate volume of fuel
        V_fuel = self.m_fuel / self.inputs.fuel.rho
        return V_fuel
    
    def V_fuel_ullage(self):
        # Calculate ullage volume in fuel tank
        V_fuel_ullage = self.inputs.fuel.tank_V - self.V_fuel()
        return V_fuel_ullage
    
    def V_cc(self):
        # Calculate combustion chamber volume
        V_cc = np.pi/4*(self.inputs.d_cc)**2*self.inputs.length_cc
        return V_cc
    
    def column_vec(self):
        # Create output column vector
        vector = np.array([
            self.m_lox,
            self.m_gox,
            self.m_oxtank_press,
            self.m_oxpresstank,
            self.m_fuel_tank_press,
            self.m_fuelpresstank,
            self.T_oxtank,
            self.T_fueltank_press,
            self.m_fuel,
            self.m_cc,
            self.M_cc,
            self.gamma_cc,
            self.T_cc
        ])
        column_vector = np.expand_dims(vector, axis=1)
        if column_vector.shape[1] == 1: # Ensure this is a column vector
            return column_vector
        else:
            return np.transpose(column_vector)

def lsv_from_column_vec(column_vec: np.ndarray, inputs: Struct):
    """ Create a LiquidStateVector object from a column vector. """
    obj                  = LiquidStateVector(inputs)
    obj.m_lox            = column_vec[0]
    obj.m_gox            = column_vec[1]
    obj.m_oxtank_press   = column_vec[2]
    obj.m_oxpresstank    = column_vec[3]
    obj.m_fueltank_press = column_vec[4]
    obj.m_fuelpresstank  = column_vec[5]
    obj.T_oxtank         = column_vec[6]
    obj.T_fueltank_press = column_vec[7]
    obj.m_fuel           = column_vec[8]
    obj.m_cc             = column_vec[9]
    obj.M_cc             = column_vec[10]
    obj.gamma_cc         = column_vec[11]
    obj.T_cc             = column_vec[12]
    obj.n2o_props        = n2o_properties(obj.T_oxtank)
    return obj


def init_liquid_state(inputs: Struct) -> Tuple[LiquidStateVec, np.ndarray]:
    """ Initializes the state vector for a liquid system. """
    # Uses the inputs to create an initial state vector for a liquid
    state_0 = LiquidStateVector(inputs),
    state_0 = init_ox_tank(state_0, inputs),
    state_0 = init_fuel_tank(state_0, inputs),
    state_0 = init_comb_chamber(state_0, inputs),
    x0 = state_0.column_vec(),
    return state_0, x0


def init_ox_tank(state: LiquidStateVector, inputs: Struct) -> LiquidStateVector:
    """ Initialize oxidizer tank properties using input Struct. """
    state.T_oxtank = inputs.ox.T_tank
    state.n2o_props = n2o_properties(state.T_oxtank)
    state.m_lox = state.n2o_props.rho_l*inputs.ox.V_l
    state.m_gox = state.n2o_props.rho_g*state.V_ox_ullage
    if inputs.ox_pressurant.active:
        state.m_oxtank_press = (inputs.ox_pressurant.set_pressure - state.p_gox)* \
            state.V_ox_ullage/(inputs.ox_pressurant.gas_props.R_specific*state.T_oxtank)
        state.m_oxpresstank = inputs.ox_pressurant.storage_init_press* \
            inputs.ox_pressurant.storage_tank_V/(inputs.ox_pressurant.gas_props.R_specific*inputs.T_amb)
    else:
        state.m_oxtank_press = 0
        state.m_oxpresstank = 0
    return state


def init_fuel_tank(state: LiquidStateVector, inputs: Struct) -> LiquidStateVector:
    """ Initialize fuel tank properties using input Struct. """
    state.T_fueltank_press = inputs.T_amb
    state.m_fuel = inputs.fuel.V_l*inputs.fuel.rho
    if inputs.fuel_pressurant.active:
        state.m_fueltank_press = inputs.fuel_pressurant.set_pressure* \
            state.V_fuel_ullage/(inputs.fuel_pressurant.gas_props.R_specific*state.T_fueltank_press)
        state.m_fuelpresstank = inputs.fuel_pressurant.storage_initial_pressure* \
            inputs.fuel_pressurant.storage_tank_V/(inputs.fuel_pressurant.gas_props.R_specific*inputs.T_amb)
    else:
        state.m_fueltank_press = 0
        state.m_fuelpresstank = 0
    return state


def init_comb_chamber(state: LiquidStateVector, inputs: Struct) -> LiquidStateVector:
    """ Initialize combustion chamber properties using input Struct. """
    M_air = 0.02897 # Mean molecular mass of air [kg/mol]
    rho_air = 1.225 # kg/m^3
    gamma_air = 1.4
    # The starting gas properties in the chamber is simply ambient air
    state.m_cc = state.V_cc*rho_air
    state.M_cc = M_air
    state.gamma_cc = gamma_air
    state.T_cc = inputs.T_amb
    return state
