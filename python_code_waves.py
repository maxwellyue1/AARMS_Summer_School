import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def seirx_model(seirx, t, parameters):
    s, e, i, r, x = seirx

    # transformed_x represents the portion of population that is willing to comply with NPIs
    transformed_x = 0.5 * np.tanh(parameters['re'] * x) + 0.5   # apply the coeff of effectiveness, then transform to [0,1]
    lambdat = parameters['beta_naught'] * (1- transformed_x)
    ds = parameters['pi'] - lambdat * s * i - parameters['mu'] * s
    de = lambdat * s * i - (parameters['kappa'] + parameters['mu']) * e
    di = parameters['kappa'] * e - (parameters['gamma'] + parameters['mu']) * i
    dr = parameters['gamma'] * i - parameters['mu'] * r
    dx = - parameters['rf'] + parameters['beta_naught'] * i

    return [ds, de, di, dr, dx]


def solve_seirx_model(re = 1, rf = 0.01): 
    '''
    re: the coefficient of effectiveness of NPIs
    rf: the level of population fatigue to NPIs
    we assume that the other parameters fixed, and S_0 = 0.995, E_0 = 0, I_0 = 0.005, R_0 = 0, X_0 = 0.8
    simulating for 1000 days, with sampling rate = 1 day
    '''



    parameters = {'pi':2e-6, 'mu':2e-6, 'beta_naught': 1.1, 'kappa': 1/6, 'gamma': 1/6, 're': re, 'rf': rf}

    s0, e0, i0, r0 = 0.995, 0, 0.005, 0
    transformed_x0 = 0.8
    x0 = np.arctanh(2 * transformed_x0 - 1) / parameters['re']  # transform back to [-inf, inf]
    print(0.5 * np.tanh(parameters['re'] * x0) + 0.5)


    start_time = 0
    end_time = 1000
    num_points = 1001
    seirx0 = [s0, e0, i0, r0, x0]  # Initial conditions
    t = np.linspace(start_time, end_time, num_points)
    solution_seirx_model = odeint(seirx_model, seirx0, t, args=(parameters, ))
    return solution_seirx_model



#t = np.linspace(0, 1000, 1001)

# plt.figure(1)
# plt.plot(t, solve_seirx_model()[:, 0], label = 'S')
# plt.plot(t, solve_seirx_model()[:, 1], label = 'E')
# plt.plot(t, solve_seirx_model()[:, 2], label = 'I')
# plt.plot(t, solve_seirx_model()[:, 3], label = 'R')
# #plt.plot(t, solve_seirx_model()[:, 4], label = 'X')
# plt.plot(t, 0.5 * np.tanh(solve_seirx_model()[:, 4]) + 0.5, label = 'transformed x')
# plt.legend()

# plt.show()
