#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brute
import ast

def get_value(input_value):
    if input_value <-55 or input_value>=32:
        return 'A1'
    elif input_value >=-55 and input_value<0:
        return 'L12'
    elif input_value >=0 and input_value<32:
        return 'L10'
    else:
        return 'bug'