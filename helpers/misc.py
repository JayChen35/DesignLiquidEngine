# Miscellaneous helper methods for DesignLiquidEngine
# Authors: Jason Chen
# Project Caelus, Aphlex 1C Engine
# Current revision: 26 February, 2021


import numpy as np
import os
import shlex
import struct
import platform
import subprocess


def print_header(string: str):
    """
    Provides an easy way to print out a header.
    :param string: The header as a string.
    """
    try:
        max_len_str = max(string.split("\n"), key=lambda x: len(x))
        size_x, size_y = get_terminal_size()
        max_length = min(len(max_len_str), size_x)
        border = "".join(["-" for _ in range(max_length)])
        print(border + "\n" + string + "\n" + border)
    except:
        print(string)


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

 
def get_terminal_size():
    """
    Gets the width and height of the current console. Works on Linux, macOS, Windows, Cygwin.
    :DISCLAIMER: This is not my code. Code was taken from the following links:
    http://stackoverflow.com/questions/566746/how-to-get-console-window-width-in-python
    https://gist.github.com/jtriley/1108174
    """
    current_os = platform.system()
    tuple_xy = None
    if current_os == 'Windows':
        tuple_xy = _get_terminal_size_windows()
        if tuple_xy is None:
            tuple_xy = _get_terminal_size_tput()
            # Needed for Window's Python in Cygwin's Xterm!
    if current_os in ['Linux', 'Darwin'] or current_os.startswith('CYGWIN'):
        tuple_xy = _get_terminal_size_linux()
    if tuple_xy is None:
        print("default")
        tuple_xy = (80, 25)
    return tuple_xy


def _get_terminal_size_windows():
    try:
        from ctypes import windll, create_string_buffer
        h = windll.kernel32.GetStdHandle(-12)
        csbi = create_string_buffer(22)
        res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)
        if res:
            (bufx, bufy, curx, cury, wattr,
             left, top, right, bottom,
             maxx, maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
            sizex = right - left + 1
            sizey = bottom - top + 1
            return sizex, sizey
    except:
        pass
 

def _get_terminal_size_tput():
    try:
        cols = int(subprocess.check_call(shlex.split('tput cols')))
        rows = int(subprocess.check_call(shlex.split('tput lines')))
        return (cols, rows)
    except:
        pass


def _get_terminal_size_linux():
    def ioctl_GWINSZ(fd):
        try:
            import fcntl
            import termios
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
            return cr
        except:
            pass
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            return None
    return int(cr[1]), int(cr[0])
