# Preliminary enigne sizing from PropSimEX2
# Authors: Jason Chen
# Project Caelus, Aphlex 1C Engine
# Current revision: 17 March, 2021

import numpy as np


def prelim_main(data: dict) -> dict:
    # Constants
    targets = dict()
    ft_to_m = 1/3.2808 # [m/ft] Conversion factor from feet to meters
    lbm_to_kg = 1/2.2046
    g_0 = 9.80665 # [m/s^2] Gravitational acceleration constant
    drag_coeff = 1.2 # [~] % more delta_v expected counter drag through flight
    dry_mass = data["mass_dry_rocket"]*lbm_to_kg
    # General parameters
    req_delta_v = np.sqrt(data["target_alt"]*ft_to_m*2*g_0) * drag_coeff
    init_mass = np.exp(req_delta_v/(data["est_isp"]*g_0)) * dry_mass
    prop_mass = init_mass-dry_mass
    impulse = req_delta_v * dry_mass # Total impulse = delta_momentum
    # Launch rail calculations
    accel = (data["min_rail_vel"]*ft_to_m)**2/(2*data["rail_length"]*ft_to_m) + g_0
    thrust = init_mass*accel
    burn_time = impulse/thrust
    # Set to data dictionary ("x" denotes "target")
    targets["x_delta_v"]   = req_delta_v
    targets["x_prop_mass"] = prop_mass
    targets["x_impulse"]   = impulse
    targets["x_thrust"]    = thrust
    targets["x_burn_time"] = burn_time
    return targets
