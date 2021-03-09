# Plotting output data from PropSimEX2
# Authors: Jason Chen
# Project Caelus, Aphlex 1C Engine
# Current revision: 05 March, 2021

import numpy as np
from scipy.io import loadmat
import matplotlib
import matplotlib.pyplot as plt


def plot_output(output_file_path: str):
    data = loadmat(output_file_path)
