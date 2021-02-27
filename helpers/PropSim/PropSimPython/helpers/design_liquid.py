# Liquid Design for PropSimPython
# Project Caelus, Aphlex 1C Engine
# Tanmay Neema, 23 February, 2021

"""
Generates a liquid engine design matching the performance
specifications and constraints provided.
   Inputs:
       - goal.max_thrust: initial thrust required
       - goal.OF: OF ratio desired
       - goal.total_impulse: total impulse desired
       - goal.min_fuel_dp: minimum pressure drop from fuel tank to
       injector as proportion of fuel tank pressure
       - goal.min_ox_dp: minimum pressure drop from oxidizer tank to
       injector as proportion of oxidizer tank pressure
       - goal.ox_to_fuel_time: ratio of liquid oxidizer flow time to liquid fuel flow time
       - design.p_tanks: initial tank pressures
       - design.ox_ullage: initial oxidizer tank volume ullage fraction
       - design.exp_ratio: design expansion ratio
   Outputs:
       - inputs: engine design parameters (see PerformanceCode.m)
"""

import numpy as np
import matplotlib.pyplot as plt
from classes import Struct
from n2o import n2o_find_T


# TODO: The following function (design_liquid) is UNTESTED
def design_liquid(initial_inputs: Struct, goal: Struct, design: Struct, output_on: bool):
    parameter_labels = ['Thrust','OF','Total Impulse','Min. Fuel Pressure Drop', 'Min. Ox. Pressure Drop','Ox. to Fuel Drain Time Ratio']
    paramater_field_names = goal.fieldnames()

    test_data = Struct()
    test_data.test_plots_on = 0 # Import tests data and plot against simulation data
    test_data.test_data_file = '' # File from which to import test data
    test_data.t_offset = 0 # Time offset of test data wrt simulation data [s]
    

    mode = Struct()
    mode.combustion_on = 1
    mode.flight_on = 0
    mode.type = 'liquid'

    max_iter = 100
    param_tol = 0.005
    
    # Enforce design parameters
    inputs = initial_inputs
    inputs.ox.T_tank = n2o_find_T(design.p_tanks)
    inputs.ox.V_tank = inputs.ox.V_l/(1-design.ox_ullage)
    inputs.fuel_pressurant.set_pressure = design.p_tanks
    inputs.exp_ratio = design.exp_ratio
