from numpy import *

# Form :
# V_plus V_minus T V_p_plus V_p_minus I_T_plus I_T_minus
results = array([[39.84e-3,39.77e-3,80.,14.24e-3,13.93e-3,858.0e-6,-58.3e-6, 80],
                 [30.37e-3,30.35e-3,80.,10.89e-3,10.66e-3,713.0e-6,-24.8e-6, 80],
                 [20.07e-3,20.06e-3,80.,7.19e-3,7.04e-3,497.9e-6,43.1e-6, 80],
                 [9.98e-3,9.97e-3,80.0,3.58e-3,3.49e-3,343.1e-6,101.4e-6, 80],

                 [9.98e-3, 9.97e-3,110.,3.94e-3,3.47e-3,800.0e-6,512.8e-6, 110],
                 [20.27e-3,20.25e-3,110.,7.71e-3,7.27e-3,813.0e-6,407.7e-6, 110],
                 [30.42e-3,30.40e-3,110.,11.40e-3,11.09e-3,850.0e-6,81.2e-6, 110],
                 [40.14e-3,40.38e-3,110.,15.14e-3,14.70e-3,1046.0e-6,40.7e-6, 110],

                 [40.38e-3,40.40e-3,31.,13.45e-3,13.17e-3,717e-6,-13.4e-6,30],
                 [30.43e-3,30.44e-3,31.,10.12e-3,9.96e-3,508.e-6,-56.8e-6,30],
                 [20.09e-3,20.09e-3,30.,6.67e-3,6.60e-3,360.2e-6,-67.5e-6,30],
                 [9.95e-3,9.94e-3,30.,3.29e-3,3.27e-3,223.6e-6,29.4e-6,30],

                 [9.95e-3,9.94e-3,50.,3.61e-3,3.09e-3,97.9e-6,-97.8e-6,50],
                 [20.08e-3,20.07e-3,50.,7.05e-3,6.62e-3,264.4e-6,-80.9e-6,50],
                 [30.51e-3,30.47e-3,50.,10.51e-3,10.17e-3,736.0e-6,-3.3e-6,50],
                 [40.21e-3,40.14e-3,50.,13.80e-3,13.47e-3,859e-6,40.3e-6,50]])

results_errors = array([[0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 80],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 80],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,0.1e-6,0.1e-6, 80],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,0.1e-6,0.1e-6, 80],

                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 110],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 110],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 110],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 110],

                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 30],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 30],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,0.1e-6,0.1e-6, 30],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,0.1e-6,0.1e-6, 30],

                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,0.1e-6,0.1e-6, 50],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,0.1e-6,0.1e-6, 50],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 50],
                        [0.01e-3,0.01e-3,1.0,0.01e-3,0.01e-3,1.e-6,0.1e-6, 50]])
