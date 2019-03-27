'''
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
fitted_pdf = norm.pdf(x_vec, loc = params[0], scale = params[1])

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.bar(x_vec, dist)
ax2.plot(x_vec, fitted_pdf)

plt.show()
'''

