# Injector parameter calculation script (like-on-like doublet impingment injector)
# Authors: Jason Chen, Ron Nachum, Jessica Chen, Tanmay Neema,
# Project Caelus, 04 March 2021

""" 
INPUTS:
    - mdot = Mass flow rate, kg/sec
    - of_ratio = Oxidizer to fuel ratio, dimensionless
    - rho_f = Fuel density, kg/m^3
    - rho_o = Oxidizer density, kg/m^3
    - P0, Chamber pressure, Pa
    - delta_p = Pressure drop across injector, % (of chamber pressure)
    - d_o = Starting diameter of oxidizer orifice, mm (1.58)
    - d_f = Starting diameter of fuel orifice, mm (1.00)
    - Cd_o = Discharge coefficient of oxidizer orifice, dimensionless (0.9)
    - Cd_f = Discharge coefficient of fuel orifice, dimensionless (0.88)
    - imp_angle = Impingement angle, degrees (60)

OUTPUTS:
    - n_o = Number of oxidizer orifices
    - d_o = Diameter of oxidizer orifice, mm
    - a_o = Area of oxidizer orifice, mm
    - L_jet_o = Oxidizer jet length, mm
    - L_poi_o = Oxidizer point of impigement distance, mm
    - d_com_o = Oxidizer orifice distance (combustor), mm
    - d_man_o = Oxidizer orifice distance (manifold), mm

    - n_f = Number of fuel orifices
    - d_f = Diameter of fuel orifice, mm
    - a_f = Area of fuel orifice, mm
    - L_jet_f = Fuel jet length, mm
    - L_poi_f = Fuel point of impigement distance, mm
    - d_com_f = Oxidizer fuel distance (combustor), mm
    - d_man_f = Oxidizer fuel distance (manifold), mm

    - L_inj = Injector plate thickness
"""

import numpy as np
import os


def injector_main(data: dict) -> dict:
    """ Calculates injector parameters. """

    mdot = data["x_mdot"]
    of_ratio = data["of_ratio"]
    mdot_o = mdot * of_ratio/(of_ratio+1)
    mdot_f = mdot * 1/(of_ratio+1)
    rho_f = data["rho_f"]
    rho_o = data["rho_o"]
    P0 = data["P0"]
    delta_p = data["delta_p"]
    d_o = data["min_d_o"]
    d_f = data["min_d_f"]
    Cd_o = data["ox"]["Cd_injector"]
    Cd_f = data["fuel"]["Cd_injector"] 
    imp_angle = data["imp_angle"]
    M = data["M_coeff"]
    jet_LD = data["jet_LD"]
    orifice_LD = data["orifice_LD"]
    
    """
    First, find the desired total injector area using an estimated Cd and target mass 
    flow rate (x_mdot), assuming constant density for nitrous oxide. Note that real nitrous
    behaves through two-phase flow (both gaseous and liquid), so increasing the injector area
    will increase the actual oxidizer mass flow rate.
    """

    if data["ox"]["injector_area"] is None or data["fuel"]["injector_area"] is None:
        # Total injector area
        A_inj_total_o = mdot_o/(Cd_o * np.sqrt(2*rho_o*P0*(delta_p/100))) 
        A_inj_total_f = mdot_f/(Cd_f * np.sqrt(2*rho_f*P0*(delta_p/100)))

    else: # Both injector areas are driving parameters
        A_inj_total_o = data["ox"]["injector_area"]
        A_inj_total_f = data["fuel"]["injector_area"]

    # Area of a single orifice [m^2]
    a_o = (np.pi * ((d_o/2) * 0.001)**2) # Attempt to use the minimum diameter orifice
    a_f = (np.pi * ((d_f/2) * 0.001)**2)
    # Number of orifices
    n_o = A_inj_total_o / a_o
    n_f = A_inj_total_f / a_f
    # Want to round down (np.floor()) so that the minimum diameter isn't exceeded
    n_o = np.floor(n_o) if np.floor(n_o) % 2 == 0 else np.floor(n_o) - 1
    n_f = np.floor(n_f) if np.floor(n_f) % 2 == 0 else np.floor(n_f) - 1
    # Diameter of a single orifice (this is different from d_o after a set n_o is chosen) [mm]
    d_o = 2 * np.sqrt(mdot_o/(Cd_o * n_o * np.pi * np.sqrt(2 * rho_o * P0 * delta_p/100))) * 1000
    d_f = 2 * np.sqrt(mdot_f/(Cd_f * n_f * np.pi * np.sqrt(2 * rho_f * P0 * delta_p/100))) * 1000
    # Length of fluid jets [mm]
    L_jet_o = jet_LD * d_o
    L_jet_f = jet_LD * d_f
    # Point of impingment [mm]
    L_poi_o = L_jet_o * np.cos(np.deg2rad(imp_angle/2))
    L_poi_f = L_jet_f * np.cos(np.deg2rad(imp_angle/2))
    # Length (thickness) of injector [mm]
    L_inj = orifice_LD * max(d_o, d_f) * np.cos(np.deg2rad(imp_angle / 2))
    # Distance between orifice (in an element pair) on combustion chamber side [mm]
    d_com_f = 2 * L_jet_f * np.sin(np.deg2rad(imp_angle / 2)) 
    d_com_o = 2 * L_jet_o * np.sin(np.deg2rad(imp_angle / 2))
    # Distance between orifice (in a pair) on injector side [mm]
    d_man_f = (d_com_o/L_poi_o) * (L_inj + L_poi_o)
    d_man_o = (d_com_f/L_poi_f) * (L_inj + L_poi_f)

    data["ox"]["injector_area"] = A_inj_total_o
    data["fuel"]["injector_area"] = A_inj_total_f
    data["x_mdot_o"] = mdot_o
    data["x_mdot_f"] = mdot_f
    data["n_o"] = n_o
    data["n_f"] = n_f
    data["d_o"] = d_o
    data["d_f"] = d_f
    data["a_o"] = a_o
    data["a_f"] = a_f
    data["L_jet_o"] = L_jet_o
    data["L_jet_f"] = L_jet_f
    data["L_poi_o"] = L_poi_o
    data["L_poi_f"] = L_poi_f
    data["L_inj"] = L_inj
    data["d_com_f"] = d_com_f
    data["d_com_o"] = d_com_o
    data["d_man_f"] = d_man_f
    data["d_man_o"] = d_man_o

    # Tank volumes = target_mdot * burn_time
    if data["ox"]["V_l"] is None: 
        data["ox"]["V_l"] = data["x_mdot_o"]*data["t_b"]
    if data["fuel"]["V_l"] is None:
        data["fuel"]["V_l"] = data["x_mdot_f"]*data["t_b"]
    
    return data


def total_cv(valve_cvs: list or np.ndarray or float) -> float:
    """ Returns the total flow coefficient (Cv) given several Cv values. """
    if type(valve_cvs) == float and valve_cvs != 0:
        return valve_cvs
    elif (type(valve_cvs) == float and valve_cvs == 0) or len(valve_cvs) == 0:
        return 0
    else:
        temp_sum = 0
        for i in range(len(valve_cvs)):
            temp_sum += 1/(valve_cvs[i]**2)
        cv = np.sqrt(temp_sum**-1)
        return cv


"""
# Mass flow rate of a single orifice
mdot_o_orifice = Cd_o * (np.pi * ((d_o/2) * 0.001)**2) * np.sqrt(2 * rho_o * P0 * delta_p/100)
mdot_f_orifice = Cd_f * (np.pi * ((d_f/2) * 0.001)**2) * np.sqrt(2 * rho_f * P0 * delta_p/100)
# Number of orifices
n_o = mdot_o / mdot_o_orifice
n_f = mdot_f / mdot_f_orifice
n_o = round(n_o) if round(n_o) % 2 == 0 else round(n_o) + 1
n_f = round(n_f) if round(n_f) % 2 == 0 else round(n_f) + 1
# Diameter of orifice
d_o = 2 * np.sqrt(mdot_o/(Cd_o * n_o * np.pi * np.sqrt(2 * rho_o * P0 * delta_p/100))) * 1000 
d_f = 2 * np.sqrt(mdot_f/(Cd_f * n_f * np.pi * np.sqrt(2 * rho_f * P0 * delta_p/100))) * 1000
# Area of orifices
a_o = np.pi * (d_o/2)**2 
a_f = np.pi * (d_f/2)**2
# Length of fluid jets
L_jet_o = jet_LD * d_o
L_jet_f = jet_LD * d_f
# Point of impingment
L_poi_o = L_jet_o * np.cos(np.deg2rad(imp_angle/2))
L_poi_f = L_jet_f * np.cos(np.deg2rad(imp_angle/2))
# Length (thickness) of injector
L_inj = orifice_LD * max(d_o, d_f) * np.cos(np.deg2rad(imp_angle / 2))
# Distance between orifice (in an element pair) on combustion chamber side
d_com_f = 2 * L_jet_f * np.sin(np.deg2rad(imp_angle / 2)) 
d_com_o = 2 * L_jet_o * np.sin(np.deg2rad(imp_angle / 2))
# Distance between orifice (in a pair) on injector side
d_man_f = (d_com_o/L_poi_o) * (L_inj + L_poi_o)
d_man_o = (d_com_f/L_poi_f) * (L_inj + L_poi_f) 
"""
