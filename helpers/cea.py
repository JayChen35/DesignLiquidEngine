# Python interface to CEA (Chemical Equilibrium with Applications)
# Authors: Tanmay Neema, Jason Chen, Daniel DeConti, Jessica Chen, Liam West
# Original interface code written by Jason Chen
# Project Caelus, 04 February 2021


import os
import shutil
import sys
import csv
import math
import numpy as np
import time
import subprocess
import threading
import pyautogui
from typing import Tuple
from helpers.misc import print_header


def cea_inp(data: dict, case_name: str) -> str:
    """ Creates CEA .inp input file given a case_name. Returns path to case directory."""
    case_dir = f"./case-files/{case_name}"
    pressure_ratio = round(data["P_0"]/data["P_3"], 3)
    # Frozen tends to overestimate engine performance, equilibrium tends to underestimate
    frozen = "frozen nfz=1" if bool(data["frozen"]) else "equilibrium" # Freezing point at throat
    # Make case folder and input file
    if os.path.exists(case_dir) and os.path.isdir(case_dir): # Clear existing folder if it exists
        shutil.rmtree(case_dir)
    os.makedirs(case_dir)

    with open(f"{case_dir}/{case_name}.inp", "w+") as f:
        line_1 = f"problem  case={case_name} o/f={data['of_ratio']},\n"
        f.write(line_1)
        line_2 = f"      rocket  {frozen}\n"
        f.write(line_2)
        line_3 = f"  p,bar={data['P_0']/1e05},\n" # Pascal to Bar conversion
        f.write(line_3)
        line_4 = f"  pi/p={pressure_ratio},\n"
        f.write(line_4)
        line_5 = "react\n"
        f.write(line_5)
        line_6 = f"  fuel=C2H5OH(L) wt=95  t,k={data['T_amb']}\n"
        f.write(line_6)
        line_7 = f"  fuel=H2O(L) wt=5  t,k={data['T_amb']}\n"
        f.write(line_7)
        line_8 = f"  oxid=N2O wt=100  t,k={data['T_amb']}\n"
        f.write(line_8)
        line_9 = "end"
        f.write(line_9)
        f.close()
    print(f"CEA input generation complete. {case_name}.inp saved to {case_dir}.")
    return case_dir


def driver_cea(case_name: str, case_dir: str):
    """ Driver for running the CEA executable. """
    cea_exe_dir = "./CEAexec/cea-exec/"
    # Move .inp file to same directory as FCEA2m.exe
    shutil.move(f"{case_dir}/{case_name}.inp", f"./CEAexec/cea-exec/{case_name}.inp")
    t1 = threading.Thread(target=type_with_delay, args=(case_name, 0.25), daemon=True)
    t1.start()
    subprocess.call([cea_exe_dir + "FCEA2m.exe"], cwd=cea_exe_dir)
    t1.join()
    # Move .inp and .out file to original directory
    shutil.move(f"./CEAexec/cea-exec/{case_name}.inp", f"{case_dir}/{case_name}.inp")
    shutil.move(f"./CEAexec/cea-exec/{case_name}.out", f"{case_dir}/{case_name}.out")
    print(f"CEA execution complete. {case_name}.out saved to {case_dir}.")

    
def type_with_delay(dir_name: str, delay: int or float):
    time.sleep(delay)
    pyautogui.typewrite(dir_name)
    pyautogui.press("enter")
   

def cea_outparse(data: dict, case_name: str, case_dir: str) -> dict:
    """ Parses the CEA .out output file to collect thermodynamics data. """
    cea_exe_dir = "./CEAexec/cea-exec/"
    csv_filename = f"{case_dir}/{case_name}.csv"
    cea_filename = f"{case_dir}/{case_name}.out"
    delimiter = "THEORETICAL ROCKET PERFORMANCE"
    with open(csv_filename, mode="w", newline="") as csv_f:
        with open(cea_filename, mode="r") as cea_f:
            cea_lines = cea_f.readlines()
            file_writer = csv.writer(csv_f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            unique_rows = set()
            row = ["O/F", "P0 (BAR)", "P1 (BAR)", "T0 (K)", "M (1/n)", "GAMMA", "CSTAR (M/SEC)", "ISP (M/SEC)"]
            dims = len(row)
            for i, line in enumerate(cea_lines): 
                if delimiter in line: # Clear the row and append to a new one
                    assert len(row) == dims # Forces the number of columns to match (a valid table)
                    if str(row) not in unique_rows: # Allows intake of multiple chamber pressures & O/F ratios
                        unique_rows.add(str(row))
                        file_writer.writerow(row)
                    row = []
                    continue
                line = line.strip()
                if "O/F=" in line:
                    temp = line.split()
                    row.append(temp[1])
                if "P, BAR" in line:
                    temp = line.split()
                    row.append(temp[2])
                    row.append(temp[3])
                    data["p_1"] = float(temp[3])*1e05 # Nozzle (critical) pressure [Pa]
                if "T, K" in line:
                    temp = line.split()
                    row.append(temp[2])
                    data["T_0"] = float(temp[2]) # Combustion chamber temperature [K]
                if "M, (1/n)" in line:
                    temp = line.split()
                    row.append(temp[2])
                    data["M"] = float(temp[2]) # Molecular mass [kg/mol]
                if "GAMMAs" in line:
                    temp = line.split()
                    row.append(temp[1])
                    data["gamma"] = float(temp[1]) # Ratio of specific heats [~]
                if "CSTAR, M/SEC" in line:
                    temp = line.split()
                    row.append(temp[2])
                    data["cstar_cea"] = float(temp[2]) # Maximum (theoretical) C* [m/s]
                if "Isp, M/SEC" in line:
                    temp = line.split()
                    row.append(temp[2])
                    data["v_2_cea"] = float(temp[2]) # Maximum (theoretical) exhaust velocity [m/s]
            file_writer.writerow(row)
    print(f"CEA output parsing complete. {case_name}.csv saved to {case_dir}.")
    return data


def cea_main(data: dict, case_name: str) -> Tuple[dict, str]:
    case_dir = cea_inp(data, case_name)
    driver_cea(case_name, case_dir)
    data = cea_outparse(data, case_name, case_dir)
    return data, case_dir
