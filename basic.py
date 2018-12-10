# -*- coding: utf-8 -*-
from numpy import *
from results import *
from thermocouple import *
from constants_and_converters import *
from mean_uncertainty import *
from matplotlib.pyplot import *

sorted_temperatures = ['30C', '50C', '80C', '110C']

T_worst_case = 2.0

temperatures = {'30C':[30.0], '50C':[50.0], '80C':[80.0], '110C':[110.0]}
temperatures_errors = {'30C':[30.0], '50C':[50.0], '80C':[80.0], '110C':[110.0]}

def select_quantity(temperature_data, temperature_data_errors, selected_quantity):
    quantity = []
    quantity_errors = []
    for i, event in enumerate(temperature_data):
        event_errors = temperature_data_errors[i]
        if i == 0:
            continue
        quantity.append(event[selected_quantity])
        quantity_errors.append(event_errors[selected_quantity])
    return array(quantity), array(quantity_errors)
    

for i, event in enumerate(results):
    batch_temperature = event[-1]
    
    event_error = results_errors[i]
    V_plus = event[0]
    s_V_plus = event_error[0]
    I_plus = V2I(V_plus)
    s_I_plus = s_V2I(s_V_plus)

    V_minus = event[1]
    s_V_minus = event_error[1]
    I_minus = V2I(V_minus)
    s_I_minus = s_V2I(s_V_minus)

    T = event[2]
    s_T = event_error[2]

    V_p_plus = event[3]
    s_V_p_plus = event_error[3]

    V_p_minus = event[4]
    s_V_p_minus = event_error[4]

    V_p, s_V_p = mean_uncertainty(array([V_p_plus, V_p_minus]), array([s_V_p_plus, s_V_p_minus]))

    I_T_plus = event[5]
    s_I_T_plus = event_error[5]
    V_T_plus = I_T_plus * R_2
    s_V_T_plus = s_I_T_plus * R_2
    T_plus, s_T_plus = thermocouple(V_T_plus, s_V_T_plus)
    s_T_plus = sqrt(s_T_plus**2 + T_worst_case**2)
    
    T_plus = T_plus
    Delta_T_plus = T_plus
    s_Delta_T_plus = s_T_plus


    I_T_minus = event[6]
    s_I_T_minus = event_error[6]
    V_T_minus = I_T_minus * R_2
    s_V_T_minus = s_I_T_minus * R_2
    T_minus, s_T_minus = thermocouple(V_T_minus, s_V_T_minus)
    s_T_minus = sqrt(s_T_minus**2 + T_worst_case**2)
    
    T_minus = T_minus
    Delta_T_minus = T_minus
    s_Delta_T_minus = s_T_minus

    difference_Delta_T = T_plus - T_minus
    s_difference_Delta_T = sqrt(s_T_plus**2 + s_T_minus**2)

    addition_Delta_T = T_plus + T_minus
    s_addition_Delta_T = sqrt(s_T_plus**2 + s_T_minus**2)

    I, s_I = mean_uncertainty(array([I_plus, I_minus]), array([s_I_plus, s_I_minus]))

    #theoretic_resistance = (rho_copper(T + abs_zero) * length_copper / F_copper + rho_constantan(T + abs_zero) * length_constantan / F_constantan)
    #s_theoretic_resistance = sqrt(s_rho_copper(s_T) **2 * length_copper ** 2 / F_copper ** 2 + s_rho_constantan(s_T) ** 2 * length_constantan ** 2 / F_constantan ** 2)
    
    #seebeck = ((V_p_plus + V_p_minus) - 2.0 * theoretic_resistance * I) / difference_Delta_T
    #s_seebeck = sqrt((s_V_plus**2 + s_V_minus**2 + 4.0 * (s_theoretic_resistance**2 * I**2 + theoretic_resistance**2 * s_I**2)) / difference_Delta_T**2 + seebeck**2 / difference_Delta_T**2 * s_difference_Delta_T)

    #V_p_replacer = ((V_p_plus + V_p_minus) - seebeck * difference_Delta_T) / 2.0
    #s_V_p_replacer = sqrt(s_V_p_plus**2 + s_V_p_minus**2 + s_seebeck**2 * difference_Delta_T**2 + seebeck**2 * s_difference_Delta_T**2)/2.0

    V_p_replacer = V_p
    s_V_p_replacer = s_V_p
   
    geometric_constant = (kappa_copper(T + abs_zero) * F_copper / length_copper + kappa_constantan(T + abs_zero) * F_constantan / length_constantan)
    s_geometric_constant = sqrt(s_kappa_copper(s_T) **2 * F_copper ** 2 / length_copper ** 2 + s_kappa_constantan(s_T) ** 2 * F_constantan ** 2 / length_constantan ** 2)

    peltier_I = geometric_constant * difference_Delta_T / (2.0 * I)
    s_peltier_I = sqrt(s_geometric_constant**2 * (difference_Delta_T / (2.0 * I))**2 + geometric_constant**2 * difference_Delta_T**2 * s_I**2 /(2.0 * I**2)**2 + geometric_constant**2 * s_difference_Delta_T**2 / (2.0 * I)**2)
    
    peltier_V = V_p_replacer / 2.0 * difference_Delta_T / addition_Delta_T
    s_peltier_V = sqrt(s_V_p_replacer**2 / 4.0 * (difference_Delta_T / addition_Delta_T)**2 + V_p_replacer**2 / 4.0 * s_difference_Delta_T**2 / addition_Delta_T**2 + V_p_replacer**2 / 4.0 * difference_Delta_T**2 / addition_Delta_T**4)

    #P_I_S = peltier_I / seebeck + abs_zero
    #s_P_I_S = sqrt(s_peltier_I**2 / seebeck**2 + peltier_I**2 / seebeck**4 * s_seebeck**2)

    #P_V_S = peltier_V / seebeck + abs_zero
    #s_P_V_S = sqrt(s_peltier_V**2 / seebeck**2 + peltier_V**2 / seebeck**4 * s_seebeck**2)
    
    #theoretic_V_p = I * theoretic_constant
    #s_theoretic_V_p = sqrt(theoretic_constant**2 * s_I**2 + s_theoretic_constant**2 * I ** 2)

    #theoretic_peltier_V = theoretic_V_p / 2.0 * difference_Delta_T / addition_Delta_T
    #s_theoretic_peltier_V = sqrt(s_theoretic_V_p**2 / 4.0 * (difference_Delta_T / addition_Delta_T)**2 + theoretic_V_p**2 / 4.0 * s_difference_Delta_T**2 / addition_Delta_T**2 + theoretic_V_p**2 / 4.0 * difference_Delta_T**2 / addition_Delta_T**4)

    #epsilon_T = V_p * difference_Delta_T / peltier_I - addition_Delta_T / 2.0
    #s_epsilon_T = sqrt(s_V_p**2 * difference_Delta_T**2 / peltier_I**2 + V_p**2 * s_difference_Delta_T**2 / peltier_I**2 + s_peltier_I**2 * V_p**2 * difference_Delta_T / peltier_I**4 + s_addition_Delta_T**2 / 2.0)
    
    for measurement in temperatures:
        if batch_temperature == int(temperatures[measurement][0]):
            temperatures[measurement].append({'I+':I_plus, 'I-':I_minus, 'I':I, 'V_p':V_p, 'DeltaT+':Delta_T_plus, 'DeltaT-':Delta_T_minus, 'DeltaT+ - DeltaT-':difference_Delta_T, 'DeltaT+ + DeltaT-':addition_Delta_T, 'peltier_I':peltier_I, 'peltier_V':peltier_V})
            temperatures_errors[measurement].append({'I+':s_I_plus, 'I-':s_I_minus, 'I':s_I, 'V_p':s_V_p, 'DeltaT+':s_Delta_T_plus, 'DeltaT-':s_Delta_T_minus, 'DeltaT+ - DeltaT-':s_difference_Delta_T, 'DeltaT+ + DeltaT-':s_addition_Delta_T, 'peltier_I':s_peltier_I, 'peltier_V':s_peltier_V})

table_file = open('results.txt', 'w')

color_marker = [['r', 'x'], ['b', '.'], ['g', 'p'], ['k', 'o']]
color_counter = 0
table_file.write('\\begin{sidewaystable}[htb]\n')
table_file.write('\\caption{Table of important quantities that were measured}\n')
table_file.write('\\centering')
table_file.write('\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\n')

for bath_temperature in sorted_temperatures:
    data = temperatures[bath_temperature]
    temperature = data[0]
    data_errors = temperatures_errors[bath_temperature]

    I_plus, s_I_plus = select_quantity(data, data_errors, 'I+')
    I_minus, s_I_minus = select_quantity(data, data_errors, 'I-')
    I, s_I = select_quantity(data, data_errors, 'I')
    V_p, s_V_p = select_quantity(data, data_errors, 'V_p')
    Delta_T_plus, s_Delta_T_plus = select_quantity(data, data_errors, 'DeltaT+')
    Delta_T_minus, s_Delta_T_minus = select_quantity(data, data_errors, 'DeltaT-')
    difference, s_difference = select_quantity(data, data_errors, 'DeltaT+ - DeltaT-')
    addition, s_addition = select_quantity(data, data_errors, 'DeltaT+ + DeltaT-')
    peltier_I, s_peltier_I = select_quantity(data, data_errors, 'peltier_I')
    peltier_V, s_peltier_V = select_quantity(data, data_errors, 'peltier_V')
    #seebeck, s_seebeck = select_quantity(data, data_errors, 'seebeck')
    #P_I_S, s_P_I_S = select_quantity(data, data_errors, 'P_I/S')
    #P_V_S, s_P_V_S = select_quantity(data, data_errors, 'P_V/S')

    sorting = argsort(I)

    I_plus = I_plus[sorting]
    s_I_plus = s_I_plus[sorting]
    I_minus = I_minus[sorting]
    s_I_minus = s_I_minus[sorting]
    I = I[sorting]
    s_I = s_I[sorting]
    V_p = V_p[sorting]
    s_V_p = s_V_p[sorting]
    Delta_T_plus = Delta_T_plus[sorting]
    s_Delta_T_plus = s_Delta_T_plus[sorting]
    Delta_T_minus = Delta_T_minus[sorting]
    s_Delta_T_minus = s_Delta_T_minus[sorting]
    difference = difference[sorting]
    s_difference = s_difference[sorting]
    addition = addition[sorting]
    s_addition = s_addition[sorting]
    peltier_I = peltier_I[sorting]
    s_peltier_I = s_peltier_I[sorting]
    peltier_V = peltier_V[sorting]
    s_peltier_V = s_peltier_V[sorting]
    #seebeck = seebeck[sorting]
    #s_seebeck = s_seebeck[sorting]
    #P_I_S = P_I_S[sorting]
    #s_P_I_S = s_P_I_S[sorting]
    #P_V_S = P_V_S[sorting]
    #s_P_V_S = s_P_V_S[sorting]

    mean_peltier_I, s_mean_peltier_I = mean_uncertainty(array(peltier_I), array(s_peltier_I))
    mean_peltier_V, s_mean_peltier_V = mean_uncertainty(array(peltier_V), array(s_peltier_V))
    #mean_seebeck, s_mean_seebeck = mean_uncertainty(array(seebeck), array(s_seebeck))

    print 'For a temperature of ' + str(temperature) + ' degrees Celsius,'
    print 'the mean values for the different quantities are : '
    print 'Peltier_I = ' + str(mean_peltier_I*1000.) + '+-' + str(s_mean_peltier_I*1000.) + 'mV'
    print 'Peltier_V = ' + str(mean_peltier_V*1000.) + '+-' + str(s_mean_peltier_V*1000.) + 'mV'
    #print 'seebeck = ' + str(mean_seebeck*1000000.) + '+-' + str(s_mean_seebeck*1000000.) + 'uV'
    print ''
    print ''

    spacing = 30
    n_columns = 10
    
    n = sorting.shape[0]
    #table_file.write('Data for ' + str(temperature) + ' degrees Celsius:\n')
    table_file.write('\\hline\n\\multicolumn{' + str(n_columns) + '}{|l|}{Results for an oil bath temperature of ' + str(temperature) + ' $^{\\circ}$C}\\\\ \\hline\\hline\n')
    
                 
    #table_file.write(' ' + ' ' * spacing * n_columns + ' ' * (n_columns-1) + ' ' + '\n')
    #table_file.write(' ' + ('=' * spacing + ' ') * (n_columns-1) + '=' * spacing + ' ' + '\n')
    table_file.write(' ' + '$I^+$'.center(spacing) + '&' + '$I^-$'.center(spacing) + '&' + '$I$'.center(spacing) + '&' + '$V_p$'.center(spacing) + '&' + '$\\Delta T^+$'.center(spacing) + '&' + '$\\Delta T^-$'.center(spacing) + '&' + '$\\Delta T^+ - \\Delta T^-$'.center(spacing) + '&' + '$\\Delta T^+ + \\Delta T^-$'.center(spacing) + '&' + '$\\Pi_{12}^{(I)}$'.center(spacing) + '&' + '$\\Pi_{12}^{(V)}$'.center(spacing) + '\\\\' + '\n')
    #table_file.write(' ' + ('=' * spacing + ' ') * (n_columns-1) + '=' * spacing + ' ' + '\n')

    for i in xrange(n):
        #table_file.write('|' + (' ' * spacing + '|') * n_columns + '\n')
        table_file.write(' ' + ('{0:.4}$\\pm${1:.4} \\si{{\\ampere}}'.format(I_plus[i], s_I_plus[i])).center(spacing) + '&' + ('{0:.4}$\\pm${1:.4} \\si{{\\ampere}}'.format(I_minus[i], s_I_minus[i])).center(spacing) + '&' + ('{0:.4}$\\pm${1:.4} \\si{{\\ampere}}'.format(I[i], s_I[i])).center(spacing) + '&' + ('{0:.4}$\\pm${1:.4} \\si{{\\milli\\volt}}'.format(V_p[i]*1000., s_V_p[i]*1000.)).center(spacing) + '&' + ('{0:.4}$\\pm${1:.4} $^{{\\circ}}$C'.format(Delta_T_plus[i], s_Delta_T_plus[i])).center(spacing) + '&' + ('{0:.4}$\pm${1:.4} $^{{\\circ}}$C'.format(Delta_T_minus[i], s_Delta_T_minus[i])).center(spacing) + '&' + ('{0:.4}$\pm${1:.4} $^{{\\circ}}$C'.format(difference[i], s_difference[i])).center(spacing) + '&' + ('{0:.4}$\pm${1:.4} $^{{\\circ}}$C'.format(addition[i], s_addition[i])).center(spacing) + '&' + ('{0:.4}$\\pm${1:.4} \\si{{\\milli\\volt}}'.format(peltier_I[i]*1000., s_peltier_I[i]*1000.)).center(spacing) + '&' + ('{0:.4}$\\pm${1:.4} \\si{{\\milli\\volt}}'.format(peltier_V[i]*1000., s_peltier_V[i]*1000.)).center(spacing) + '\\\\' + '\n')
        #table_file.write('|' + ('_' * spacing + '|') * n_columns + '\n')

    #table_file.write('\n\n\n')
    table_file.write('\\hline\n')

    figure(0) # Plot of DeltaT+
    errorbar(I, Delta_T_plus, xerr = s_I, yerr = s_Delta_T_plus, color = color_marker[color_counter][0], marker = color_marker[color_counter][1], linestyle = ':', label = r'$T_b$ = ' + bath_temperature[:-1] + r'$^{\circ}$C')
    figure(1) # Plot of DeltaT-
    errorbar(I, Delta_T_minus, xerr = s_I, yerr = s_Delta_T_minus, color = color_marker[color_counter][0], marker = color_marker[color_counter][1], linestyle = ':', label = r'$T_b$ = ' + bath_temperature[:-1] + r'$^{\circ}$C')
    figure(2) # Plot of DeltaT+ - DeltaT-
    errorbar(I, difference, xerr = s_I, yerr = s_difference, color = color_marker[color_counter][0], marker = color_marker[color_counter][1], linestyle = ':', label = r'$T_b$ = ' + bath_temperature[:-1] + r'$^{\circ}$C')
    figure(3) # Plot of DeltaT+ + DeltaT+
    errorbar(I, addition, xerr = s_I, yerr = s_addition, color = color_marker[color_counter][0], marker = color_marker[color_counter][1], linestyle = ':', label = r'$T_b$ = ' + bath_temperature[:-1] + r'$^{\circ}$C')
    figure(4) # Plot of peltier_I
    errorbar(I, peltier_I*1000., xerr = s_I, yerr = s_peltier_I*1000., color = color_marker[color_counter][0], marker = color_marker[color_counter][1], linestyle = ':', label = r'$T_b$ = ' + bath_temperature[:-1] + r'$^{\circ}$C')
    figure(5) # Plot of peltier_V
    errorbar(I, peltier_V*1000., xerr = s_I, yerr = s_peltier_V*1000., color = color_marker[color_counter][0], marker = color_marker[color_counter][1], linestyle = ':', label = r'$T_b$ = ' + bath_temperature[:-1] + r'$^{\circ}$C')

    color_counter += 1
table_file.write('\\end{tabular}\n')
table_file.write('\\end{sidewaystable}')
table_file.close()

figure(0)
gca().tick_params(bottom = True, top = True, left = True, right = True)
#grid(True)
legend(loc = 'best')
title(r"Evolution of $\Delta T^+$ for different current intensities $I$" + "\n" + r"and oil bath temperatures $T_b$", fontsize = 18)
xlabel(r'Current intensity $I$ in [A]', fontsize = 20)
ylabel(r'$\Delta T^+$ in [K]', fontsize = 20)
savefig('DeltaT_plus.pdf')

figure(1)
gca().tick_params(bottom = True, top = True, left = True, right = True)
#grid(True)
legend(loc = 'best')
title(r"Evolution of $\Delta T^-$ for different current intensities $I$" + "\n" + r"and oil bath temperatures $T_b$", fontsize = 18)
xlabel(r'Current intensity $I$ in [A]', fontsize = 20)
ylabel(r'$\Delta T^-$ in [K]', fontsize = 20)
ylim(-2.0, 6.0)
savefig('DeltaT_minus.pdf')

figure(2)
gca().tick_params(bottom = True, top = True, left = True, right = True)
#grid(True)
legend(loc = 'best')
title(r"Evolution of $\Delta T^+ - \Delta T^-$ for different current intensities $I$" + "\n" + r"and oil bath temperatures $T_b$", fontsize = 18)
xlabel(r'Current intensity $I$ in [A]', fontsize = 20)
ylabel(r'$\Delta T^+ - \Delta T^-$ in [K]', fontsize = 20)
savefig('difference.pdf')

figure(3)
gca().tick_params(bottom = True, top = True, left = True, right = True)
#grid(True)
legend(loc = 'best')
title(r"Evolution of $\Delta T^+ + \Delta T^-$ for different current intensities $I$" + "\n" + r"and oil bath temperatures $T_b$", fontsize = 18)
xlabel(r'Current intensity $I$ in [A]', fontsize = 20)
ylabel(r'$\Delta T^+ + \Delta T^-$ in [K]', fontsize = 20)
savefig('addition.pdf')

figure(4)
gca().tick_params(bottom = True, top = True, left = True, right = True)
#grid(True)
legend(loc = 'best')
title(r"Peltier coefficient $\Pi_{12}^I$ for different current intensities $I$" + "\n" + r"and oil bath temperatures $T_b$", fontsize = 18)
xlabel(r'Current intensity $I$ in [A]', fontsize = 20)
ylabel(r'$\Pi_{12}^I$ in [mV]', fontsize = 20)
savefig('peltier_I.pdf')

figure(5)
gca().tick_params(bottom = True, top = True, left = True, right = True)
#grid(True)
legend(loc = 'best')
title(r"Peltier coefficient $\Pi_{12}^V$ for different current intensities $I$" + "\n" + r"and oil bath temperatures $T_b$", fontsize = 18)
xlabel(r'Current intensity $I$ in [A]', fontsize = 20)
ylabel(r'$\Pi_{12}^V$ in [mV]', fontsize = 20)
ylim(0.0, 10.0)
savefig('peltier_V.pdf')

#show()
    
