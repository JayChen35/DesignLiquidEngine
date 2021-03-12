# DesignLiquidEngine - Written in Python 3.6+
# Designs a liquid bipropellant rocket engine given input parameters.
# Project Caelus, Aphlex 1C Engine
# Current revision: 26 February, 2021
# Authors: Jason Chen, Tanmay Neema, Ron Nachum, Anya Mischel, Jessica Chen

import numpy as np
import os
import ast
import yaml
import sys
import time
import datetime
import platform
from typing import Tuple
from helpers.misc import get_exit_pressure, print_header
from helpers.cea import cea_main
from helpers.nozzle import nozzle_main
from helpers.injector import injector_main
from helpers.propsim import propsim_main, check_valid_input
from helpers.output import compile_outputs


def main(cmd_args: list) -> Tuple[dict, float, str]:
    # Ingest user parameters and store into data{} dictionary.
    contains_cmd_args = len(cmd_args) > 1
    if len(cmd_args) == 2 and (cmd_args[1] == "--default" or cmd_args[1] == "-d"):
        temp_path, temp_name = "", ""
    elif len(cmd_args) == 3 and (cmd_args[1] == "-c" or cmd_args[1] == "-n"):
        temp_path = "" if cmd_args[1] == "-n" else cmd_args[2] # Inputting a name = default config
        temp_name = "" if cmd_args[1] == "-c" else cmd_args[2] # Inputting a config = default name
    elif len(cmd_args) == 5 and (cmd_args[1] == "-c") and (cmd_args[3] == "-n"):
        temp_path, temp_name = cmd_args[2], cmd_args[4]
    elif len(cmd_args) == 5 and (cmd_args[1] == "-n") and (cmd_args[3] == "-c"):
        temp_path, temp_name = cmd_args[4], cmd_args[2]
    elif not contains_cmd_args:
        pass
    else:
        print_header("Invalid command line arguments. Please see README.md for correct formatting.")
        sys.exit(0)
    terminal = False
    while not terminal:
        if not contains_cmd_args:
            print_header("Enter configuration (parameter) file path, including the `.yaml` extension.\
                \nAlternatively, press [ENTER] to use default file (./config.yaml), or [0] to exit.")
            temp_path = input("Please enter a file path or command: ")
        try:
            if temp_path == "":
                temp_path = "./config.yaml"
            elif temp_path == "0":
                sys.exit(0)
            with open(temp_path, "r") as f:
                data = yaml.safe_load(f)
            assert type(data) == dict
            terminal = True
            if not check_valid_input(data, {np.float64, float, int, bool, type(None), str}):
                print_header("Invalid configuration file (some data types not allowed). Try again.")
                sys.exit(0)
        except (IOError, ValueError, AssertionError):
            print_header("Invalid file path. Please ensure the file path is typed correctly.")
            if contains_cmd_args: sys.exit(0)
    # Get ambient pressure given the input altitude
    data["P3"] = get_exit_pressure(data["altitude"])
    # Add precision to data type
    data["ox"]["injector_area"] = np.float64(data["ox"]["injector_area"])
    data["fuel"]["injector_area"] = np.float64(data["fuel"]["injector_area"])
    # Get desired name of output/case files from user
    if not contains_cmd_args:
        temp_name = input("Enter the desired name of output/case files (press [ENTER] to use default): ")
    if temp_name == "":
        case_name = datetime.datetime.now().strftime("%d-%b-%Y_%H%M%S") + "EST"
    else: # Make the inputted case name a valid case name
        case_name = "".join([x if(x.isalnum() or x in "._- ") else "_" for x in temp_name])
    # Change current working directory to ./helpers
    os.chdir("./helpers")
    start_time = time.perf_counter() # Start runtime timer after inputs are recieved
    # Run Chemical Equilibrium with Applications (CEA) interface code
    data, case_dir = cea_main(data, case_name)
    # Run nozzle design code (assuming quasi-1D isentropic flow)
    data = nozzle_main(data)
    # Run injector design code (assuming 2-on-1 doublet impinging injector)
    data = injector_main(data)
    # Run PropSim (liquid engine performance simulator)
    output_file_path = propsim_main(data, case_dir)
    # Compile output data into a summary
    design_dict, _ = compile_outputs(data, output_file_path)
    # Print calculation outputs
    for key, val in design_dict.items():
        print(f"{key} = {val}")
    return data, start_time, case_dir, output_file_path


# TODO: Add graphing of pressures post-MPV in MATLAB
# TODO: Design around impulse not t_b, implement save_data_on, organize PropSimInput.json
# TODO: Implement mass gradient descent optimizer
# TODO: Add more instructions/documentation to README.md
if __name__ == "__main__":
    cmd_args = sys.argv # Accepts command line arguements
    data, start_time, case_dir, output_file_path = main(cmd_args)
    duration = round(time.perf_counter()-start_time, 4)
    case_dir = case_dir[1:].replace("/", "\\") if platform.system() == "Windows" else case_dir[1:]
    output_path = os.getcwd() + case_dir
    print_header(f"DesignLiquidEngine completed execution in {duration} seconds.\nOutputs saved to {output_path}.")
