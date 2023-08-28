import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


parameters = {'pi':2e-6, 'mu':2e-6, 'beta_naught': 1.1, 'kappa': 1/6, 'gamma': 1/6, 'k': 1
 }

s0, e0, i0, r0 = 0.995, 0, 0.005, 0
x0 = 0.8
start_time = 0
end_time = 1000
num_points = 1001

seirx0 = [s0, e0, i0, r0, x0]  # Initial conditions
t = np.linspace(start_time, end_time, num_points)

def seirx_model(seirx, t, parameters):
    s, e, i, r, x = seirx

    lambdat = parameters['beta_naught'] * (1- parameters['k'] * (0.5 * np.tanh(x) + 0.5))
    ds = - lambdat * s * i
    de = lambdat * s * i - parameters['kappa'] * e
    di = parameters['kappa'] * e - parameters['gamma'] * i
    dr = parameters['gamma'] * i
    dx = - 0.01 + parameters['beta_naught'] * i


    return [ds, de, di, dr, dx]

solution_seirx_model = odeint(seirx_model, seirx0, t, args=(parameters, ))




plt.figure(1)
plt.plot(t, solution_seirx_model[:, 0], label = 'S')
plt.plot(t, solution_seirx_model[:, 1], label = 'E')
plt.plot(t, solution_seirx_model[:, 2], label = 'I')
plt.plot(t, solution_seirx_model[:, 3], label = 'R')
plt.plot(t, 0.5 * np.tanh(solution_seirx_model[:, 4]) + 0.5, label = 'transformed x')
plt.legend()

plt.show()
