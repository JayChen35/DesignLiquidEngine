# DesignLiquidEngine - Written in Python 3.6+
# Designs a liquid bipropellant rocket engine given input parameters.
# Project Caelus, Aphlex 1C Engine
# Current revision: 26 February, 2021
# Authors: Jason Chen, Tanmay Neema, Ron Nachum, Anya Mischel, Jessica Chen, Cyril Sharma


import numpy as np
import os
import yaml
import sys
import time
import datetime
from helpers.misc import get_exit_pressure, print_header
from helpers.cea import cea_main
from helpers.nozzle import nozzle_main
from helpers.injector import injector_main


def main() -> dict:
    # Ingest user parameters and store into data{} dictionary.
    terminal = False
    while not terminal:
        print("Enter configuration (parameter) file path, including .yaml extension (e.g. \"config.yaml\").")
        print("Alternatively, press [ENTER] to use default file (./config.yaml), or [0] to exit.")
        temp_path = input("File path: ")
        try:
            if temp_path == "":
                temp_path = "./config.yaml"
            elif temp_path == "0":
                sys.exit(0)
            with open(temp_path, "r") as f:
                data = yaml.safe_load(f)
                assert type(data) == dict
                terminal = True
        except IOError or ValueError or AssertionError:
            print_header("Invalid file path. Please ensure the file path is typed correctly.")
    input_path = temp_path
    # Get ambient pressure given the input altitude
    data["P3"] = get_exit_pressure(data["altitude"])
    # Get desired name of output/case files from user
    temp_name = input("Enter the desired name of output/case files (press [ENTER] to use default): ")
    if temp_name == "":
        case_name = datetime.datetime.now().strftime("%d-%b-%Y_%H%M%S") + "EST"
    else:
        case_name = "".join([x if(x.isalnum() or x in "._- ") else "_" for x in temp_name])
    # Change current working directory to ./helpers
    os.chdir("./helpers")
    # Run Chemical Equilibrium with Applications (CEA) interface code
    data = cea_main(data, case_name)
    # Run nozzle design code (assuming quasi-1D isentropic flow)
    data = nozzle_main(data)
    # Run injector design code (assuming 2-on-1 doublet impinging injector)
    data = injector_main(data)
    # Print calculation outputs
    for key, val in data.items():
        print(f"{key} = {val}")
    return data


def plot_output(data: dict):
    pass


def save_output(data: dict):
    pass


if __name__ == "__main__":
    start = time.perf_counter()
    data = main()
    if data["plot_on"]:
        plot_output(data)
    if data["save_data_on"]:
        save_output(data)
    duration = round(time.perf_counter()-start, 4)
    print_header(f"DesignLiquidEngine completed exectution in {duration} seconds.")
