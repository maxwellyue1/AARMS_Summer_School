import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


omega_const = 100
stringencies = [100]
parameters = {'pi':2e-6, 'mu':2e-6, 'beta_naught': 1.1, 'kappa': 1/6, 'gamma': 0.125, 'delta': 0.001, 'K': 0.4 }

s0, e0, i0, r0 = 0.995, 0, 0.005, 0
x0 = 0.8
start_time = 0
end_time = 600
num_points = 601

#x = 1

seirx0 = [s0, e0, i0, r0, x0]  # Initial conditions
t = np.linspace(start_time, end_time, num_points)

solutions_diff_stringencies = {}

for stringency in stringencies: 

    def seirx_model(seirx, t, parameters):
        s, e, i, r, x = seirx
        
        lambdat = parameters['beta_naught'] * i / (1 + stringency * i * x)
        ds = parameters['pi'] - (lambdat + parameters['mu']) * s
        de = lambdat * s - (parameters['kappa'] + parameters['mu']) * e
        di = parameters['kappa'] * e - (parameters['gamma'] + parameters['delta'] + parameters['mu']) * i
        dr = parameters['gamma'] * i - parameters['mu'] * r
        dx = parameters['K'] * x * (1 - x) * (-1 + omega_const * i)
        

        return [ds, de, di, dr, dx]
    
    solution = odeint(seirx_model, seirx0, t, args=(parameters, ))
    
    solutions_diff_stringencies[f'alpha={stringency}'] = solution

for _ in solutions_diff_stringencies.keys(): 
     plt.plot(t, solutions_diff_stringencies[_][:, 2], label = _)
     plt.plot(t, parameters['beta_naught'] * (1-solutions_diff_stringencies[_][:, 4]), label = _)

# Save the plot to a file
#plt.legend()
plt.savefig('figure1.png')

# Display the plot
plt.show()