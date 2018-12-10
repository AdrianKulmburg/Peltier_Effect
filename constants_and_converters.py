from numpy import *

abs_zero = 273.15

R_1 = 1200.0
R_2 = 389e-3

rho_copper = lambda x : (1.68 + (x - 293.0) * (2.06 - 1.68) / (353.0 - 293.0))*1e-8
s_rho_copper = lambda s_x : s_x * abs((2.06 - 1.68)/(353.0 - 293.0)) * 1e-8
kappa_copper = lambda x : 401.0 + (x - 250.0) * (391.0 - 401.0)/(400.0 - 250.0)
s_kappa_copper = lambda s_x : s_x * abs((391.0 - 401.0)/(400.0 - 250.0))
length_copper = 50e-3
diameter_copper = 2.03e-3
F_copper = pi * (diameter_copper / 2.0) ** 2

rho_constantan = lambda x : 44.0*1e-8
s_rho_constantan = lambda s_x : 0.0
kappa_constantan = lambda x : 21.9 + (x - 275.0) * (26.6 - 21.9)/(400.0 - 275.0)
s_kappa_constantan = lambda s_x : s_x * abs((26.6 - 21.9)/(400.0 - 275.0))
length_constantan = 50e-3
diameter_constantan = 7.04e-3
F_constantan = pi * (diameter_constantan / 2.0) ** 2

# In order to convert from V to I
V2I = lambda v : v / 50.0e-3 * 20.0
s_V2I = lambda s_v : s_v / 50.0e-3 * 20.0
