import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
from scipy.constants import Boltzmann as kb

N = 1000
steps = 100
field = 2*steps+1



'''Random Walk in potential'''
particles = np.zeros(N)
k = 0.5e-20
h = 1
T=300
beta = 1/(kb*T)

def P_min(x, V):
    rel_prob = np.exp(-beta*(V(x-h)-V(x+h)))
    Pp = 1/(1 + rel_prob)
    Pm = 1 - Pp
    return Pm

def rand_walk_pot(p_pos, V):
    for i in range(steps):
        steps_vec = np.random.rand(N)
        for x in p_pos:
            steps_vec[steps_vec >= P_min(x, V)] = 1
            steps_vec[steps_vec != 1] = -1
        p_pos += steps_vec
        params = norm.fit(p_pos)
    return p_pos, params

# Define potentials
V0 = lambda x: 0
V1 = lambda x: k*x

def V2(x):
    if -3*h < x and x < 3*h:
        return k
    else:
        return 0

p, params = rand_walk_pot(particles, V1)
x_vec = np.arange(min(p), max(p))
fitted_pdf = norm.pdf(x_vec, loc = params[0], scale = params[1])
print(params)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.hist(p)
ax2.plot(x_vec, fitted_pdf, color="r")
plt.show()
