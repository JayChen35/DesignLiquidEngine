# DesignLiquidEngine - Written in Python 3.6+
# Designs a liquid bipropellant rocket engine given input parameters.
# Project Caelus, Aphlex 1C Engine
# Current revision: 26 February, 2021
# Authors: Jason Chen, Tanmay Neema, Ron Nachum, Anya Mischel, Jessica Chen, Cyril Sharma


import numpy as np
import os
from helpers.cea import cea_main, get_exit_pressure
from helpers.nozzle import nozzle_main
from helpers.injector import injector_main
import yaml


def take_all_inputs():
    print("Input variables using a .txt file. Enter floats only. Ensure file is within folder of .py file.")
    print("Please list them in the order they are listed at the top of the file, i.e. \n Thrust=___ \n Chamber Pressure=___ \n . \n . \n .")
    
    file_ext = input("Enter file path with .txt: ")
    file_name = file_ext.split(".")[0]
    data = dict()
    try:
        with open(file_ext) as f:
            for line in f:
                equal_index = line.find("=")
                name = line[:equal_index].strip()
                value = float(line[equal_index+1:].strip())
                data[name] = value
    except IOError:
        print("Please ensure input.txt is named correctly and in the correct directory.")
    except ValueError:
        print("Please ensure inputs are entered as floats with no other text in the file")
    data["P3"] = get_exit_pressure(data["altitude"])
    ceagui_name = file_name + "_ceagui"
    return data, ceagui_name


def print_header(string, key=lambda: 30):
    """
    Provides an easy way to print out a header.
    :param string: The header as a string.
    :param key: Length of the message.
    """
    print("\n")
    print((len(string) % 2)*'-' + '{:-^{width}}'.format(string, width=key()))


if __name__ == "__main__":
    py_dir = os.path.dirname(__file__)
    # os.chdir(py_dir)
    os.chdir("./helpers")
    temp = take_all_inputs()
    data = temp[0]
    ceagui_name = temp[1]
    cea_main(data, ceagui_name)
    nozzle_main(data)
    injector_main(data)
    for key in data:
        print(f"{key} = {data[key]}")
