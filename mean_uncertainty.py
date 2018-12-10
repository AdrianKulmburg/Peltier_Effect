from numpy import *

def mean_uncertainty(values, values_errors):
    return sum(values / values_errors ** 2) / sum(1.0 / values_errors ** 2), sqrt(1.0 / sum(1.0 / values_errors ** 2))
