# Supersonic Nozzle Design Calculation Program
# Created by Jason Chen, Project Caelus
# Thomas Jefferson High School for Science and Technology
# Developed on April 22nd, 2018

""" 

This program generates important parameters in a conical supersonic converging-diverging
nozzle from certain user inputs. Below is a full list of this program's capabilities and
key thermodynamic assumptions. Many input parameters can be obtained from NASA's CEA program
(Chemical Equilibrium with Applications), which is available freely available to the
public for download here (https://www.grc.nasa.gov/www/CEAWeb/).


-------------------------------------------------------------------------------------
ALL CALCULATIONS ARE DONE IN STANDARD IMPERIAL (SI) UNITS, WITH TEMPERATURE IN KELVIN.
ANGLES IN DEGREES.
-------------------------------------------------------------------------------------


CALCULATION ASSUMPTIONS:
    1) Flow maintains equilibrium at the local pressure and enthalpy level (reversible).
    2) No heat transfer across gas-enclosure walls; closed system (adiabatic).
    3) By definition, the flow is isentropic due to assumptions 1 and 2.
    4) Flow composition is considered frozen (tends to underestimate performance by 1-4%).
    5) Working fluid obeys the perfect gas law.
    6) All species of the working fluid in the flow is gaseous.
    7) Wall friction is negligible and all boundary layer effects are ignored.
    8) There are no shock waves or other discontinuities within the nozzle flow.
    9) Propellant mass flow rate is constant and steady.
    10) Expansion of the working fluid is uniform and steady; no significant turbulence.
    11) All exhaust gases travel with a velocity parallel to the nozzle axis.
    12) Gas velocity, pressure, temperature, & density are all uniform across any section
        perpendicular to the nozzle axis.
    13) Nozzle dimensions assume an optimum divergence half-angle of 15 degrees, and an
        optimum convergence half-angle of 45 degrees.
    14) Contraction ratio (cross-sectional area of the chamber/area of the throat) of 8.

INPUTS:
    - F = Desired thrust, N
    - P0 = Chamber pressure, Pa
    - ALT = Altitude, m
    - O/F = Oxidizer to fuel ratio, dimensionless
    - T0 = Combustion chamber temperature, K
    - M = Molecular mass of the gas, kg/mol
    - k = Ratio of specific heat capacities, cp/cv, dimensionless
    - L* = Characteristic chamber length, m

OUTPUTS:
    - Isp = Specific impulse at altitude, sec
    - Tt = Throat temperature, K
    - v2 = Effective exhaust velocity, m/sec
    - mdot = Mass flow rate, kg/sec
    - mdot_oxidizer = Mass flow rate of the oxidizer, kg/sec
    - mdot_fuel = Mass flow rate of the fuel, kg/sec
    - PR = Pressure ratio, dimensionless
    - ER = Expansion ratio, dimensionless
    - Te = Exit temperature, K
    - Mnum = Exit Mach number, dimensionless
    - At = Area of the throat, m^2
    - Ae = Area of the exit, m^2
    - Rt = Radius of the throat, m
    - Re = Radius of the exit, m
    - Rc = Radius of the chamber, m
    - Lc = Length of the chamber, m
    - Ldn = Length of the diverging nozzle, m
    - Lcn = Length of the converging nozzle, m

"""

# Imports
import numpy as np
import os



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
    data["P3"] = get_exit_pressure(data["altitude"])
    return data

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




def calculate(data):
    """
    Attempts to calculate and print values.
    """
    F = data["thrust"]
    P0 = data["P0"]
    ALT = data["altitude"]
    P3 = data["P3"]
    OF = data["o/f"]
    T0 = data["T0"]
    M = data["M"]
    k = data["GAMMA"]
    Lstar = data['Lstar']

    try:  # Attempt to calculate values
        R = (8314.3 / M)
        PR = (P3 / P0)
        AR = (((k + 1) / 2) ** (1 / (k - 1))) * ((P3 / P0) ** (1 / k)) * (
            np.sqrt(((k + 1) / (k - 1)) * (1 - ((P3 / P0) ** ((k - 1) / k)))))
        ER = 1 / AR
        Tt = (2 * T0) / (k + 1)
        v2 = np.sqrt((2 * k / (k - 1)) * ((R) * T0) * (1 - ((P3 / P0) ** ((k - 1) / k))))
        mdot = F / v2
        mdot_fuel = (mdot / (OF + 1))
        mdot_oxidizer = (mdot / (OF + 1)) * OF
        Isp = F / (mdot * 9.80655)
        Te = T0 / ((P0 / P3) ** ((k - 1) / k))
        Mnum = (v2 / (np.sqrt(k * (R) * (Te))))
        At = ((mdot) * (np.sqrt((k * R * T0)))) / (k * P0 * (np.sqrt(((2 / (k + 1)) ** ((k + 1) / (k - 1))))))
        Ae = ER * At
        Rt = np.sqrt(At / np.pi)
        Re = np.sqrt(Ae / np.pi)
        Ac = At * 8
        Rc = np.sqrt((Ac) / np.pi)
        Lc = ((At) * Lstar) / (Ac)
        Ldn = ((Re) - (Rt)) / (np.tan(np.deg2rad(15)))
        Lcn = ((Rc) - (Rt)) / (np.tan(np.deg2rad(45)))

        data["R"] = R
        data["PR"] = PR
        data["AR"] = AR
        data["ER"] = ER
        data["Tt"] = Tt
        data["v2"] = v2
        data["mdot"] = mdot
        data["mdot_fuel"] = mdot_fuel
        data["mdot_oxidizer"] = mdot_oxidizer
        data["Isp"] = Isp
        data["Te"] = Te
        data["Mnum"] = Mnum
        data["At"] = At
        data["Ae"] = Ae
        data["Rt"] = Rt
        data["Re"] = Re
        data["Ac"] = Ac
        data["Rc"] = Rc
        data["Lc"] = Lc
        data["Ldn"] = Ldn
        data["Lcn"] = Lcn

    except (ValueError, ZeroDivisionError):  # Exception thrown
        print("\n", "Error while attempting to solve. Please enter a valid value for every parameter.")

def print_outputs(data):
    for key in data:
        print(f"{key} = {data[key]}")


def nozzle_main(data):
    calculate(data)



if __name__ == "__main__":
    py_dir = os.path.dirname(__file__)
    os.chdir(py_dir)
    data = input_variables()
    calculate(data)
    print_outputs(data)
