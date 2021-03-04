# Miscellaneous helper methods for DesignLiquidEngine
# Authors: Jason Chen
# Project Caelus, Aphlex 1C Engine
# Current revision: 26 February, 2021


import numpy as np


def print_header(string: str):
    """
    Provides an easy way to print out a header.
    :param string: The header as a string.
    """
    border = "".join(["-" for _ in range(len(string))])
    print(border + "\n" + string + "\n" + border)


def print_seperator(string: str, key=lambda: 30):
    """
    Provides an easy way to print out a header.
    :param string: The header as a string.
    :param key: Length of the message.
    """
    print("\n")
    print((len(string) % 2)*'-' + '{:-^{width}}'.format(string, width=key()))


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
        P3 = (22.65 * np.exp(1.73 - 0.000157 * h))*1000
    else: # Troposphere
        T = 15.04 - 0.00649 * h
        P3 = (101.29 * ((T + 273.1) / 288.08) ** (5.256))*1000
    return P3


class Struct():
    """
    Generic data class to store named variables. Can be initialized recursively using an input dictionary
    to work with nested dictionary structures. Used to easily call variables with the "." operator.
    """
    def __init__(self, *args, **kwargs):
        super().__init__()
        if len(args) == 1 and type(args[0]) == dict: # Create from a dictionary
            for key, value in args[0].items():
                if isinstance(value, dict): # If the value is a dictionary (more nests to deal with)
                    setattr(self, key, Struct(value)) # Recursively create another Struct instance
                else:
                    setattr(self, key, value) # If the value is regular data, set the key, value pair
    
    def __getitem__(self, key): # Defines the "." operator to access a value in the internal dictionary.
        return self.__dict__[key]

    def __repr__(self): # __repr__ allows printing of this class to be more readible.
        return "{%s}" % str(", ".join("%s : %s" % (k, repr(v)) for (k, v) in self.__dict__.items()))

    def get_dict(self):
        return self.__dict__
