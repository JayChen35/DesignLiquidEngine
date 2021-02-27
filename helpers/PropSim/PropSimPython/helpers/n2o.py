# Nitrous-related helper methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# Liam West, Anya Mischel, & Jason Chen, 10 February, 2021


import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from classes import Struct


def n2o_properties(temp: int or float) -> Struct:
    """
    Calculates an array of properties of nitrous oxide given a temperature in K.
    WARNING: if temperature input is outside of -90 to 30 C, properties will
    be generated for boundary (-90 C or 30 C).
    """
    properties = dict() # Creates properties (output array) as a dictionary for which data from text file can be entered
    properties["Pvap"]       = None # Nitrous vapor pressure
    properties["rho_l"]      = None # Liquid density of nitrous
    properties["rho_g"]      = None # Gas denisty of nitrous
    properties["deltaH_vap"] = None # Vaportization enthalpy
    properties["cp_l"]       = None # Specific heat of liquid nitrous
    properties["cv_l"]       = None # Specific volume of liquid nitrous
    properties["cp_g"]       = None # Specific heat of gaseous nitrous
    properties["cv_g"]       = None # Specific volume of gaseous nitrous 
    properties["h_l"]        = None # Enthalpy of liquid nitrous
    properties["h_g"]        = None # Enthalpy of gaseous nitrous
    properties["s_l"]        = None # Specific entropy of liquid nitrous
    properties["s_g"]        = None # Specific entropy of gaseous nitrous
    properties["mu_l"]       = None # Dynamic viscosity of nitrous liquid
    properties["mu_g"]       = None # Dynamic viscosity of nitrous gas
    properties["e_l"]        = None # Specific internal energy of liquid
    properties["e_g"]        = None # Specific internal energy of gas
        
    R_u = 8.3144621 # Universal gas constant [J/mol*K]
    M_n2o = 0.044013 # Molecular mass of nitrous oxide [kg/mol]
    R_n2o_0 = R_u/M_n2o # Specific gas constant of nitrous oxide [J/kg*K]

    # Range-check temperature
    if temp < (-90 + 273.15): # If temperature is less that bottom bound
        temp = -90 + 273.150001  # Temperature is bottom bound for interpolate
    elif temp > (30 + 273.150001): # If temperature greater than top bound, 
        temp = 30 + 273.150001 # Temperature equal to top bound for interpolate

    Tcrit = 309.57 # K
    Pcrit = 7251 # kPa
    rhocrit = 452 # kg/m^3
    # Possibly add critical compressibility factor "Z"
    Tr = temp/Tcrit

    # Calculate vapor pressure, valid -90 to 36C
    b1 = -6.71893
    b2 = 1.3596
    b3 = -1.3779
    b4 = -4.051

    properties["Pvap"] = np.exp((1/Tr)*(b1*(1-Tr) + b2*(1-Tr)**(3/2) + b3*(1-Tr)**(5/2) + b4*(1-Tr)**5))*Pcrit
    properties["Pvap"] = properties["Pvap"]*1000

    # Calculate Density of Liquid, valid -90C to 36C
    b1 = 1.72328
    b2 = -0.83950
    b3 = 0.51060
    b4 = -0.10412

    properties["rho_l"] = np.exp(b1*(1-Tr)**(1/3) + b2*(1-Tr)**(2/3) + b3*(1-Tr) + b4*(1-Tr)**(4/3))*rhocrit

    # Calculate Density of Gas, valid -90C to 36C
    b1 = -1.00900
    b2 = -6.28792
    b3 = 7.50332
    b4 = -7.90463
    b5 = 0.629427
    Trinv = 1./Tr

    properties["rho_g"] = np.exp(b1*(Trinv-1)**(1/3) + b2*(Trinv-1)**(2/3) + b3*(Trinv-1) + b4*(Trinv-1)**(4/3) + b5*(Trinv-1)**(5/3))*rhocrit

    # Calculate dynamic viscosity of saturated liquid, valid from -90C to 30C
    b1 = 1.6089
    b2 = 2.0439
    b3 = 5.24
    b4 = 0.0293423
    theta = (Tcrit-b3)/(temp-b3)

    properties["mu_l"] = b4*np.exp(b1*(theta-1)**(1/3) + b2*(theta-1)**(4/3))

    # Calculate dynamic viscosity of saturated vapor, valid from -90C to 30C
    b1 = 3.3281
    b2 = -1.18237
    b3 = -0.055155
    Trinv = 1./Tr

    properties["mu_g"] = np.exp(b1 + b2*(Trinv-1)**(1/3) + b3*(Trinv-1)**(4/3))
    NIST_data = dict() # NIST_data is an array that stores variables regarding nitrous

    # Read each line in the N2O_Properties.cgi.txt document and enter each line into the dictionary, separated by tabs
    with open("./data/N2O_Properties.cgi.txt", "r") as reader:
        for x, line in enumerate(reader):
            if x ==0:
                temp_list = line.split("\t")
                arr = [list() for x in range(len(temp_list))]
                continue
            temp_list = line.split("\t") 
            for i, val in enumerate(temp_list):
                arr[i].append(val)
        NIST_data["T"] =     arr[0]
        NIST_data["h_liq"] = arr[5]
        NIST_data["h_gas"] = arr[17]
        NIST_data["e_liq"] = arr[4]
        NIST_data["e_gas"] = arr[16]
        NIST_data["cv_l"]  = arr[7]
        NIST_data["cv_g"]  = arr[19]
        NIST_data["s_l"]   = arr[6]
        NIST_data["s_g"]   = arr[18]
        NIST_data["cp_liq"]= arr[8]
        NIST_data["cp_gas"]= arr[20]
    NIST_data = {key: np.array(NIST_data[key], dtype="float64") for key in NIST_data}

    # Gas Specific Enthalpy
    properties["h_l"] = np.interp(temp, NIST_data["T"], NIST_data["h_liq"])*1000 # J/kg,   
    properties["h_g"] = np.interp(temp, NIST_data["T"], NIST_data["h_gas"])*1000 # J/kg   
    properties["e_l"] = np.interp(temp, NIST_data["T"], NIST_data["e_liq"])*1000 # J/kg   
    properties["e_g"] = np.interp(temp, NIST_data["T"], NIST_data["e_gas"])*1000 # J/kg   
    properties["deltaH_vap"] = properties["h_g"]-properties["h_l"]
    properties["deltaE_vap"] = properties["e_g"]-properties["e_l"]

    # Specific Heat at Constant Volume
    properties["cv_l"] = np.interp(temp, NIST_data["T"], NIST_data["cv_l"])*1000
    properties["cv_g"] = np.interp(temp, NIST_data["T"], NIST_data["cv_g"])*1000

    # Specific Heat at Constant Pressure
    properties["cp_l"] = np.interp(temp, NIST_data["T"], NIST_data["cp_liq"])*1000
    properties["cp_g"] = np.interp(temp, NIST_data["T"], NIST_data["cp_gas"])*1000

    # Specific Entropy
    properties["s_l"] = np.interp(temp, NIST_data["T"], NIST_data["s_l"])*1000
    properties["s_g"] = np.interp(temp, NIST_data["T"], NIST_data["s_g"])*1000

    # Convert Properties to Standard Units
    properties["mu_l"] = properties["mu_l"]*10**-3 # mN*s/(m^2) -> N*s/m^2
    properties["mu_g"] = properties["mu_g"]*10**-6 # uN*s/(m^ 2)-> N*s/m^2

    props = Struct(properties) # Converts the output properties into standardized Struct type
    return props


def get_n2o_Pvap(temp: int or float):
    """A faster method to get the nitrous vapor pressure based on input temperature."""

    # Range-check temperature
    if temp < (-90 + 273.15): # If temperature is less that bottom bound
        temp = -90 + 273.150001  # Temperature is bottom bound for interpolate
    elif temp > (30 + 273.150001): # If temperature greater than top bound, 
        temp = 30 + 273.150001 # Temperature equal to top bound for interpolate

    Tcrit = 309.57 # K
    Pcrit = 7251 # kPa
    rhocrit = 452 # kg/m^3
    Tr = temp/Tcrit

    # Calculate vapor pressure, valid -90 to 36C
    b1 = -6.71893
    b2 = 1.3596
    b3 = -1.3779
    b4 = -4.051

    Pvap = 1000*(np.exp((1/Tr)*(b1*(1-Tr) + b2*(1-Tr)**(3/2) + b3*(1-Tr)**(5/2) + b4*(1-Tr)**5))*Pcrit)
    return Pvap


def n2o_find_T(p_vap: int or float) -> float:
    """ Estimate temperature that gives the specified vapor pressure. """
    T_min = -90 + 273.15
    T_max = 30 + 273.15
    num_pts = 1000
    T_i = np.linspace(T_min, T_max, num_pts)
    Pvap = np.zeros(T_i.shape)

    for x in range(num_pts):
        Pvap[x] = get_n2o_Pvap(T_i[x])

    z = InterpolatedUnivariateSpline(Pvap, T_i, k=1) # Best-fit line, Pvap along x, temp along y
    return z(p_vap)


def two_phase_n2o_flow(t_1: int or float):    
    """
    Calculates the critical mass flow rate using the homogeneous equilibrium 
    model. The mass flow rate may not exceed this, regardless of the back pressure.

    Inputs:
    T_1: upstream liquid temperature, K

    Outputs:
       crit_flow:
           crit_flow.p_up_norm: vector of (upstream pressure / upstream vapor 
               pressure) for G_crit calculations
           crit_flow.G_crit: critical oxidizer mass flux at each p_1, kg/m^2*s
           crit_flow.p_down_norm_crit: critical normalized back pressure 
               (downstream pressure / upstream pressure) at each p_1, Pa
       mass_flux:
           mass_flux.p_up_norm: array of (upstream pressure / upstream vapor 
               pressure)
           mass_flux.p_down_norm: array of (back pressure / upstream pressure)
           mass_flux.G: array of mass flux at corresponding P_1_norm, P_2_norm

    NOTE: States are denoted by subscript as follows:
        1: upstream (stagnation)
        2: downstream
    """

    # Numerical options
    n = 200 # Density of calculation points along each dimension

    # Saturation properties at upstream temperature
    n2o_prop_1_sat = n2o_properties(t_1)

    # Create array of p_1, p_2
    p_1_min = 1.0
    p_1_max = 3.0
    p_1 = np.linspace(p_1_min, p_1_max, n)
    p_2 = np.linspace(0, 1, n)

    P_1_norm, P_2_norm = np.meshgrid(p_1, p_2)

    # Convert to pressures (from pressure ratio)
    P_1 = np.multiply(P_1_norm, n2o_prop_1_sat.Pvap)
    P_2 = np.multiply(P_2_norm, P_1)
    
    # Upstream liquid enthalpy
    # Calculate enthalpy as extension from saturated state
    enthalpy_1 = n2o_prop_1_sat.h_l + (P_1 - n2o_prop_1_sat.Pvap)/n2o_prop_1_sat.rho_l

    # Calculate downstream temperature
    T_2 = n2o_find_T(P_2)
    # n2o_prop_2 = # [n2o_properties(j) for j in i for i in T_2]
    a = T_2.shape[0] 
    b = T_2.shape[1]
    n2o_prop_2 = [[0 for b in range(b)] for i in range(a)]

    for x in range(a):
        for y in range(b):
            n2o_prop_2[x][y] = n2o_properties(T_2[x, y])

    # Calculate gas density
    rho_2_g = np.array([[j.rho_g for j in i] for i in n2o_prop_2])
    rho_2_l = n2o_prop_1_sat.rho_l * np.ones(np.array(P_2).shape)

    enthalpy_2_g = np.array([[j.h_g for j in i] for i in n2o_prop_2])
    enthalpy_2_l = np.array([[j.h_l for j in i] for i in n2o_prop_2])

    entropy_1_l = n2o_prop_1_sat.s_l * np.ones(np.array(P_2).shape)
    entropy_2_l = np.array([[j.s_l for j in i] for i in n2o_prop_2])
    entropy_2_g = np.array([[j.s_g for j in i] for i in n2o_prop_2])

    # entropy difference: liquid downstream - liquid upstream
    entropy_ld_lu = entropy_2_l - entropy_1_l
    # entropy difference: liquid downstream - gas downstream
    entropy_ld_gd = entropy_2_l - entropy_2_g

    # calculate mass fraction of gas to conserve entropy
    # massfrac: mass fraction of vapor
    massfrac = np.divide(entropy_ld_lu, entropy_ld_gd)

    a = massfrac.shape[0] 
    b = massfrac.shape[1]
    for x in range(a):
        for y in range(b):
            if massfrac[x, y] < 0:
                rho_2_l = n2o_prop_1_sat.rho_l
                enthalpy_2_l = n2o_prop_1_sat.h_l + P_2 - n2o_prop_1_sat.Pvap / n2o_prop_1_sat.rho_l
                massfrac[x, y] = 0

    # downstream inverse density
    rho_2_inv = np.multiply(massfrac, (1./rho_2_g) + np.multiply((1-massfrac), (1./rho_2_l)))

    # downstream equivalent density
    rho_2_equiv = np.ones(rho_2_inv.shape) / rho_2_inv

    enthalpy_2 = massfrac * enthalpy_2_g + (np.ones(massfrac.shape) - massfrac) * enthalpy_2_l

    # Homogeneous Equilibrium Model
    G = rho_2_equiv * np.sqrt(2 * (enthalpy_1 - enthalpy_2))

    # G_crit for each 
    G_crit = np.amax(G, axis=0) 
    i_crit = np.argmax(G, axis=0)

    P_2_crit = np.zeros((1, P_2.shape[1]))
    for ii in range(0, P_2.shape[1]):
        P_2_crit[ii] = P_2[i_crit[ii], ii]

    # Create downstream pressure vs. oxidizer mass flux profile
    G_out = np.matlib.repmat(G_crit, [P_2.shape[0], 1])
    # P_2_crit expanded to size of P_2
    P_2_crit_exp = np.matlib.repmat(P_2_crit, [P_2.shape[0], 1])
    G_out = G if P_2 > P_2_crit_exp else None

    # Calculate incompressible mass flux
    G_inc = np.sqrt(2*(P_1-P_2)*n2o_prop_1_sat.rho_l)

    # Package for output
    crit_flow = {}
    mass_flux = {}
    crit_flow["p_up_norm"] = p_1
    crit_flow["G_crit"] = G_crit
    crit_flow["p_down_norm_crit"] = np.divide(P_2_crit, (p_1*n2o_prop_1_sat.Pvap))
    mass_flux["p_up_norm"] = P_1_norm
    mass_flux["p_down_norm"] = P_2_norm
    mass_flux["G"] = G_out

    return crit_flow, mass_flux

two_phase_n2o_flow(300)
