# Compiling output data from PropSimEX2
# Authors: Jason Chen
# Project Caelus, Aphlex 1C Engine
# Current revision: 05 March, 2021

import numpy as np
import json
import scipy.signal
import matplotlib
import matplotlib.pyplot as plt
from helpers.misc import loadmat, Struct


def compile_outputs(data: dict, output_file_path: str) -> tuple:
    g_0 = 9.80665 # [m/s^2] Gravitational acceleration constant
    design = dict()
    record = loadmat(output_file_path)["record"] # Custom loadmat to solve bad formatting
    max_thrust, m_dot_max, p_cc_max = get_max_thrust(record)
    impulse = record["impulse"]
    isp = record["Isp"]/g_0
    # Mass of propellants
    Mox_initial = record["m_ox"].flatten()[0]
    Mox = record["m_ox"].flatten()[0] - record["m_ox"].flatten()[-1]
    Vox = (Mox/data["rho_o"])*1e03 # [L]
    Mfuel_initial = record["m_fuel"].flatten()[0]
    Mfuel = record["m_fuel"].flatten()[0] - record["m_fuel"].flatten()[-1]
    Vfuel = (Mfuel/data["rho_f"])*1e03 # [L]
    of_ratio = Mox/Mfuel
    avg_mdot_ox = np.mean(record["m_dot_ox"].flatten()[1:]) # First element is None
    avg_mdot_fuel = np.mean(record["m_dot_fuel"].flatten()[1:]) # First element is None
    # Propellant tank parameters
    start_p_oxtank = record["p_oxtank"].flatten()[0]
    end_p_oxtank = record["p_oxtank"].flatten()[-1]
    start_p_fueltank = record["p_fueltank"].flatten()[0]
    end_p_fueltank = record["p_fueltank"].flatten()[-1]
    # Other parameters
    exit_mach = get_filtered_max(record, record["M_e"].flatten())
    local_c = np.sqrt(1.4*287.058*data["T_amb"]) # Local speed of sound, 1.4 gamma, 287.058 J/kg*K
    burn_time = record["time"].flatten()[-1]
    wet_mass = data["mass_dry_rocket"] + Mfuel_initial + Mox_initial
    dry_mass = data["mass_dry_rocket"] + (Mfuel_initial-Mfuel) + (Mox_initial-Mox)
    delta_v = isp*g_0*np.log(wet_mass/dry_mass)
    # delta_v = exit_mach*local_c*np.log(wet_mass/dry_mass)
    ideal_alt = 0.5*(delta_v**2)/g_0 # Simply using delta_K = delta_U; no air resistance
    # Add vital data to output dictionary
    # TODO: ADD MORE CRITICAL VALUES HERE; add values from data_dict to here
    design["max_thrust"]       = max_thrust
    design["m_dot_max"]        = m_dot_max
    design["p_cc_max"]         = p_cc_max
    design["impulse"]          = impulse
    design["isp"]              = isp
    design["Mox_initial"]      = Mox_initial
    design["Mox_used"]         = Mox
    design["Vox_used"]         = Vox
    design["Mfuel_initial"]    = Mfuel_initial
    design["Mfuel_used"]       = Mfuel
    design["Vfuel_used"]       = Vfuel
    design["of_ratio"]         = of_ratio
    design["avg_mdot_ox"]      = avg_mdot_ox
    design["avg_mdot_fuel"]    = avg_mdot_fuel
    design["start_p_oxtank"]   = start_p_oxtank
    design["end_p_oxtank"]     = end_p_oxtank
    design["start_p_fueltank"] = start_p_fueltank
    design["end_p_fueltank"]   = end_p_fueltank
    design["A_inj_ox_eff"]     = data["ox"]["injector_area"] # Only the injector area (orifice sizes)
    design["A_inj_fuel_eff"]   = data["fuel"]["injector_area"]
    design["A_inj_o_only"]     = data["ox"]["A_inj_o_only"] # Effective injector area (including Cv)
    design["A_inj_f_only"]     = data["ox"]["A_inj_f_only"]
    design["exit_mach"]        = exit_mach
    design["burn_time"]        = burn_time
    design["ideal_delta_v"]    = delta_v
    design["ideal_alt"]        = ideal_alt
    # Cut off extraneous significant figures
    for key, val in design.items():
        design[key] = np.round(val, 4) if val > 0.1 else round(val, 5-int(np.floor(np.log10(abs(val))))-1)
    # Export to JSON file
    json_obj = json.dumps(design, indent=4, separators=(",", ": "))
    prefix = output_file_path[:output_file_path.rfind("/")]
    json_path = prefix + "/FinalDesignSummary.json"
    with open(json_path, "w+") as f:
        f.write(json_obj)
    return design, json_path


def get_max_thrust(record: dict) -> tuple:
    """ Filters out local spikes and finds max thrust, mdot, and chamber pressure. """
    dt_filter = 0.1
    dn_thrust_filter = np.ceil(dt_filter/np.mean(np.diff(record["time"])))
    a = np.array([1])
    b = (1/dn_thrust_filter * np.ones((int(dn_thrust_filter), 1))).flatten()
    filtered_thrust = scipy.signal.lfilter(b, a, record["F_thrust"])
    filtered_m_dot = scipy.signal.lfilter(b, a, record["m_dot_prop"])
    filtered_p_cc = scipy.signal.lfilter(b, a, record["p_cc"])
    max_ind = np.argmax(filtered_thrust)
    max_thrust = filtered_thrust[max_ind]
    m_dot_max = filtered_m_dot[max_ind]
    p_cc_max = filtered_p_cc[max_ind]
    return max_thrust, m_dot_max, p_cc_max


def get_filtered_max(record: dict, vals: list or np.ndarray) -> float:
    dt_filter = 0.1
    dn_filter = np.ceil(dt_filter/np.mean(np.diff(record["time"])))
    a = np.array([1])
    b = (1/dn_filter * np.ones((int(dn_filter), 1))).flatten()
    filtered_val = scipy.signal.lfilter(b, a, vals)
    max_ind = np.argmax(filtered_val)
    return filtered_val[max_ind]


if __name__ == "__main__":
    compile_outputs(dict(), "./case-files/test_case_2/PropSimOutput.mat")
