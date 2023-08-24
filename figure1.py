import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

stringencies = [0, 10, 50, 100, 200]
parameters = {'beta_naught': 1.1, 'kappa': 1/6, 'gamma': 0.125, 'delta': 0.001, 'K': 0.4 }

s0, e0, i0, r0 = 0.8, 0, 0.2, 0
start_time = 0
end_time = 300
num_points = 301

x = 4

seir0 = [s0, e0, i0, r0]  # Initial conditions
t = np.linspace(start_time, end_time, num_points)

def lambda_fct(i, x, parameters, stringency): 
        lambdat = parameters['beta_naught'] * i / (1 + stringency * i * x)
        return lambdat

solutions_diff_stringencies = {}

for stringency in stringencies: 

    def seir_model(seir, t, parameters):
        s, e, i, r = seir
        lambdat = lambda_fct(i, x, parameters)
        ds = - lambdat * s
        de = lambdat * s - parameters['kappa'] * e
        di = parameters['kappa'] * e - (parameters['gamma'] + parameters['delta']) * i
        dr = parameters['gamma'] * i

        return [ds, de, di, dr]
    
    solution = odeint(seir_model, seir0, t, args=(parameters, ))
    
    solutions_diff_stringencies[f'alpha={stringency}'] = solution

for _ in solutions_diff_stringencies.keys(): 
     plt.plot(t, solutions_diff_stringencies[_][:, 2])

# Save the plot to a file
plt.savefig('plot_example.png')

# Display the plot
plt.show()