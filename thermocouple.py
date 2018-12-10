from numpy import *
from general_fit import general_fit

voltages_table = array([-0.392, -0.353, -0.314, -0.275, -0.236, -0.197, -0.157, -0.118,
                  -0.079, -0.039, 0.000, 0.000, 0.039, 0.079, 0.119, 0.158, 0.198,
                  0.238, 0.277, 0.317, 0.357, 0.397, 0.397, 0.437, 0.477, 0.517, 0.557, 0.597, 0.637, 0.677, 0.718, 0.758, 0.798]) * 1e-3
temperatures_table = array([-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 0, 1, 2, 3,
                      4, 5, 6, 7, 8, 9, 10, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                      19, 20]) + 0.0

def linear_fit(parameters, data):
    a = parameters[0]
    b = parameters[1]

    return a * data + b

def thermocouple(voltage, voltage_error):
    odr_parameter_ideal, odr_parameter_error, odr_p_value, odr_SSR, odr_covariance = general_fit(voltages_table, temperatures_table, linear_fit, [3.9e-5, 0], return_covariance = True)
    a = odr_parameter_ideal[0]
    b = odr_parameter_ideal[1]
    s_a = odr_parameter_error[0]
    s_b = odr_parameter_error[1]
    x = linspace(min(voltages_table), max(voltages_table), 1000)
    y = a * x + b
    #print a, s_a
    #print b, s_b

    temperature = a * voltage + b
    temperature_error = sqrt(s_a ** 2 * voltage ** 2 + a ** 2 * voltage_error ** 2 + s_b ** 2 + 2 * voltage * odr_covariance[1, 0])

    return temperature, temperature_error
