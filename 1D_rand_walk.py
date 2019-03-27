import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm

N = 1000
steps = 100
field = 2*steps+1
x_vec = np.linspace(-steps, steps, field)

# This uses some weird matrix technic, cause I want to try to remove the inner for-loop
particles = np.zeros([N, field])
particles.T[steps] = np.ones(N)

for i in range(steps):
    steps_vec = np.around(np.random.rand(N))
    steps_vec[steps_vec == 0] = -1
    for j in range(N):
        particles[j] = np.roll(particles[j], int(steps_vec[j]))

dist = np.sum(particles, 0)
params = norm.fit(dist)
print(params)
fitted_pdf = norm.pdf(x_vec, loc = params[0], scale = params[1])

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.bar(x_vec, dist)
ax2.plot(x_vec, fitted_pdf)

plt.show()


'''Random Walk in potential'''
from scipy.constants import Boltzmann as kb

particles = np.zeros(N)
k = 0.5e-20
h = 1
T=300
beta = 1/(kb*T)
print(beta*k)
V = lambda x: k*x

def P_min(x):
    rel_prob = np.exp(-beta*(V(x-h)-V(x+h)))
    Pp = 1/(1 + rel_prob)
    Pm = 1 - Pp
    return Pm


for i in range(steps):
    steps_vec = np.random.rand(N)
    for x in particles:
        steps_vec[steps_vec >= P_min(x)] = 1
        steps_vec[steps_vec != 1] = -1
    particles += steps_vec

plt.hist(np.sort(particles), x_vec)
plt.show()
