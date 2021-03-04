# Injector parameter calculation script (2-on-1 impingment injector)
# Authors: Ron Nachum, Jason Chen, Jessica Chen, Tanmay Neema, 
# Project Caelus, 04 March 2021

""" 
-------------------------------------------------------------------------------------
ALL CALCULATIONS ARE DONE IN STANDARD IMPERIAL (SI) UNITS, WITH TEMPERATURE IN KELVIN.
ANGLES IN DEGREES.
-------------------------------------------------------------------------------------

INPUTS:
    - mdot = Mass flow rate, kg/sec
    - of_ratio = Oxidizer to fuel ratio, dimensionless
    - rho_f = Fuel density, kg/m^3
    - rho_o = Oxidizer density, kg/m^3
    - P0, Chamber pressure, Pa
    - delta_p = Pressure drop across injector, % (of chamber pressure)
    - d_o = Starting diameter of oxidizer orifice, mm (1.58)
    - Cd_o = Discharge coefficient of oxidizer orifice, dimensionless (0.9)
    - d_f = Starting diameter of fuel orifice, mm (1.00)
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
    """ Attempts to calculate and print values. """

    mdot = data["mdot"]
    of_ratio = data["of_ratio"]
    mdot_o = mdot * of_ratio/(of_ratio+1)
    mdot_f = mdot * 1/(of_ratio+1)
    rho_f = data["rho_f"]
    rho_o = data["rho_o"]
    P0 = data["P0"]
    delta_p = data["delta_p"]
    og_d_o = data["d_o"]
    Cd_o = data["Cd_o"]
    og_d_f = data["d_f"]
    Cd_f = data["Cd_f"]
    imp_angle = data["imp_angle"]
    M = data["M_coeff"]
    jet_LD = data["jet_LD"]
    orifice_LD = data["orifice_LD"]

    diam_ratio = np.sqrt(M * (((rho_o/rho_f)*(mdot_o/mdot_f)**2)**0.7))
    
    mdot_o_orifice = Cd_o * (np.pi * ((og_d_o/2) * 0.001)**2) * np.sqrt(2 * rho_o * P0 * 0.25)
    mdot_f_orifice = Cd_f * (np.pi * ((og_d_f / 2) * 0.001)**2) * np.sqrt(2 * rho_f * P0 * 0.25)
    
    n_o = mdot_o / mdot_o_orifice
    n_f = mdot_f / mdot_f_orifice
    
    n_o = round(n_o) if round(n_o) % 2 == 0 else round(n_o) + 1
    n_f = round(n_f) if round(n_o) % 2 == 0 else round(n_f) + 1

    d_o = 2 * np.sqrt(mdot_o/(Cd_o * n_o * np.pi * np.sqrt(2 * rho_o * P0 * delta_p/100))) * 1000
    d_f = 2 * np.sqrt(mdot_f/(Cd_f * n_f * np.pi * np.sqrt(2 * rho_f * P0 * delta_p/100))) * 1000

    a_o = np.pi * (d_o/2)**2
    a_f = np.pi * (d_f/2)**2

    L_jet_o = jet_LD * d_o
    L_jet_f = jet_LD * d_f

    L_poi_o = L_jet_o * np.cos(np.deg2rad(imp_angle/2))
    L_poi_f = L_jet_f * np.cos(np.deg2rad(imp_angle/2))
    
    L_inj = orifice_LD * max(d_o, d_f) * np.cos(np.deg2rad(imp_angle / 2))

    d_com_f = 2 * L_jet_f * np.sin(np.deg2rad(imp_angle / 2))
    d_com_o = 2 * L_jet_o * np.sin(np.deg2rad(imp_angle / 2))
    d_man_f = (d_com_o/L_poi_o) * (L_inj + L_poi_o)
    d_man_o = (d_com_f/L_poi_f) * (L_inj + L_poi_f)

    data["diam_ratio"] = diam_ratio
    data["mdot_o_orifice"] = mdot_o_orifice
    data["mdot_f_orifice"] = mdot_f_orifice
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
    data["d_man_o"] = d_man_o
    data["d_man_f"] = d_man_f
    data["d_man_o"] = d_man_o
    
    return data
