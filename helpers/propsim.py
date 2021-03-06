# Driver for the compiled MATLAB PropSim program
# Authors: Jason Chen
# Project Caelus, Aphlex 1C Engine
# Current revision: 05 March, 2021


import json
import sys
import shutil
import os
import numpy as np
from helpers.misc import print_header


def propsim_main(data: dict, case_dir: str):
    """ Compiles parameters and calls propsim.exe or propsim.sh. """
    # Check that the data dictionary is valid and all nulls are filled
    data["options"]["output_on"] = True if data["plot_on"] else False
    if not check_valid_input(data):
        error_string = "ERROR: The configuration (.yaml) file received invalid inputs, or there was an internal"+\
            " calculation error. Check that all nulls in the YAML are in valid places, and try running main.py"+\
            " with the default `config.yaml` file to check basic capability. See README.md for more details."
        print_header(error_string)
        sys.exit(0)
    data["case_name"] = case_name = case_dir.split("/")[-1]
    # Convert to JSON and export to the case directory. This will later be read using MATLAB.
    json_obj = json.dumps(data, indent=4, separators=(',', ': '))
    json_path = case_dir + f"/{case_name}.json"
    json_temp_path = f"./PropSim/{case_name}.json"
    with open(json_path, "w+") as f:
        f.write(json_obj)
    # Copy the JSON file to the same directory as PropSim
    shutil.copyfile(json_path, json_temp_path)
    # Call the compiled MATLAB executable to run PropSim
    pass
    # Delete the temporary JSON file
    # os.remove(json_temp_path)


def check_valid_input(data, allowable_types: set = {np.float64, float, int, bool}) -> bool:
    """ Recursively checks the validity of the current data dictionary. """
    # All values should be ints, floats, arrays, or booleans.
    if type(data) in allowable_types:
        return True
    elif type(data) == list:
        valid = True
        for item in data:
            valid = check_valid_input(item, allowable_types) and valid
        return valid
    elif type(data) == dict:
        valid = True
        for item in data.values():
            valid = check_valid_input(item, allowable_types) and valid
        return valid
    else:
        return False
