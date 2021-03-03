
""" 
-------------------------------------------------------------------------------------
ALL CALCULATIONS ARE DONE IN STANDARD IMPERIAL (SI) UNITS, WITH TEMPERATURE IN KELVIN.
ANGLES IN DEGREES.
-------------------------------------------------------------------------------------

INPUTS:
    - mdot = Mass flow rate, kg/sec
    - OF = Oxidizer to fuel ratio, dimensionless
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

# Imports
import numpy as np
import os



def get_exit_pressure(h: int or float):
    """
    Sourced from NASA Glenn Research Center's Earth Atmosphere Model.
    Computes and sets the ambient pressure, P3, based on inputted altitude (meters).
    P3 has units in pascals. Note: The intermediate temperature calculations use Celsius.
    :param h: Altitude, in meters.
    :return P3: Ambient pressure at altitude, in pascals.
    """
    if (h >= 25000):  # Upper Stratosphere
        T = -131.21 + 0.00299 * h
        P3 = (2.488 * ((T + 273.1) / 216.6) ** (-11.388))*1000
    elif (11000 < h < 25000):  # Lower Stratosphere
        T = -56.46
        P3 = (22.65 * math.exp(1.73 - 0.000157 * h))*1000
    else: # Troposphere
        T = 15.04 - 0.00649 * h
        P3 = (101.29 * ((T + 273.1) / 288.08) ** (5.256))*1000
    return P3

def input_variables():
    print("Input variables using a .txt file. Enter floats only. Ensure file is within folder of .py file.")
    print("Please list them in the order they are listed at the top of the file, i.e. \n Thrust=___ \n Chamber Pressure=___ \n . \n . \n .")
    
    file_ext = input("Enter file path with .txt: ")
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
    print(data)
    P3 = get_exit_pressure(data["altitude"])
    data["P3"] = P3
    return data




def calculate(data):
    """
    Attempts to calculate and print values.
    """
    try:  # Attempt to calculate values
        mdot = data["mdot"]
        OF = data["o/f"]
        mdot_o = mdot * OF/(OF+1)
        mdot_f = mdot * 1/(OF+1)
        rho_f = data["rho_f"]
        rho_o = data["rho_o"]
        P0 = data["P0"]
        delta_p = data["delta_p"]
        og_d_o = data["og_d_o"]
        Cd_o = data["Cd_o"]
        og_d_f = data["og_d_f"]
        Cd_f = data["Cd_f"]
        imp_angle = data["imp_angle"]
        M = data["M"]
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
    except (ValueError, ZeroDivisionError):  # Exception thrown
        print("\n", "Error while attempting to solve. Please enter a valid value for every parameter.")
    
    return data


def print_out():
    for key in data:
        print(f"{key} = {data[key]}")

def injector_main(data):
    calculate(data)


if __name__ == "__main__":
    py_dir = os.path.dirname(__file__)
    os.chdir(py_dir)
    data = input_variables()
    calculate(data)
    print_outputs(data)