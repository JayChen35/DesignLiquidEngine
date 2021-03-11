# Driver for the compiled MATLAB PropSim program
# Authors: Jason Chen
# Project Caelus, Aphlex 1C Engine
# Current revision: 05 March, 2021

import json
import sys
import shutil
import os
import numpy as np
import subprocess
import platform
from helpers.misc import print_header


def propsim_main(data: dict, case_dir: str) -> str:
    """ Compiles parameters and calls propsim.exe or propsim.sh. """
    # Check that the data dictionary is valid and all nulls are filled
    data["options"]["output_on"] = True if data["plot_on"] else False
    if not check_valid_input(data, {np.float64, float, int, bool, str}):
        error_string = "ERROR: The configuration (.yaml) file received invalid inputs, or there was an internal"+\
            " calculation error. Check that all nulls in the YAML are in valid places, and try running main.py"+\
            " with the default `config.yaml` file to check basic capability. See README.md for more details."
        print_header(error_string)
        sys.exit(0)
    data["case_name"] = case_name = case_dir.split("/")[-1]
    # Check for valid input/output file names
    input_file_name = "".join([x if(x.isalnum() or x in "._- ") else "_" for x in data["propsim_io"]["input_file"]])
    data["propsim_io"]["input_file"] = input_file_name
    output_file_name = "".join([x if(x.isalnum() or x in "._- ") else "_" for x in data["propsim_io"]["output_file"]])
    data["propsim_io"]["output_file"] = output_file_name
    # Convert to JSON and export to the case directory. This will later be read using MATLAB.
    json_obj = json.dumps(data, indent=4, separators=(",", ": "))
    json_path = case_dir + f"/{input_file_name}.json"
    json_temp_path = f"./{input_file_name}.json"
    with open(json_path, "w+") as f:
        f.write(json_obj)
    # Copy the JSON file to the same directory as PropSim
    shutil.copyfile(json_path, json_temp_path)
    # Call the compiled MATLAB executable to run PropSim
    if platform.system() == "Windows":
        subprocess.call(f"./PropSimIntegrated.exe {json_temp_path}")
    elif platform.system() in ["Linux", "Darwin", "Java"]:
        subprocess.call(f"./run_PropSimIntegrated.sh $MCR_PATH/v99 {json_temp_path}")
    # Delete the temporary JSON file and move output .mat file from PropSim to case-files
    os.remove(json_temp_path)
    output_file_path = case_dir + f"/{output_file_name}.mat"
    shutil.move(f"./{output_file_name}.mat", output_file_path)
    info = f"PropSimIntegrated.exe successfully executed. Output .mat saved to {output_file_path}."
    if data["save_data_on"]:
        output_plot_path = f"./{output_file_name}Plot.fig"
        if os.path.exists(output_plot_path):
            new_plot_path = case_dir + f"/{output_file_name}Plot.fig"
            shutil.move(output_plot_path, new_plot_path)
            info += f" Output MATLAB figure saved to {new_plot_path}."
    print_header(info)
    return output_file_path


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
