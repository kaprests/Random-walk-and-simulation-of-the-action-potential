import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm

N = 1000
steps = 100
field = 2*steps
x_vec = np.linspace(-steps, steps, field)

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
