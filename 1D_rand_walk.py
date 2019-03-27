import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
from scipy.constants import Boltzmann as kb, elementary_charge as elemc

steps = 100
particles = np.zeros(N)
h = 1
T0 = 273+37 #K 
T=T0
beta_k = 10
beta = 1/(kb*T)
k = beta_k/beta

def P_min(x, V):
    rel_prob = np.exp(-beta*(V(x-h)-V(x+h)))
    Pp = 1/(1 + rel_prob)
    Pm = 1 - Pp
    return Pm


def rand_walk(p_pos, V, N=1000, L=0):
    for i in range(steps):
        steps_vec = np.random.rand(N)
        for x in p_pos:
            steps_vec[steps_vec >= P_min(x, V)] = 1
            steps_vec[steps_vec != 1] = -1
        p_pos += steps_vec
        if L != 0:
            steps_vec[steps_vec == L+1] = L-1
            steps_vec[steps_vec == -L-1] = -L+1
        params = norm.fit(p_pos)
    return p_pos, params


def plot_distribution(p_pos, V):
    p, params = rand_walk(p_pos, V)
    x_vec = np.arange(min(p), max(p))
    fitted_pdf = norm.pdf(x_vec, loc = params[0], scale = params[1])
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.hist(p, bins=50)
    ax2.plot(x_vec, fitted_pdf, color="r")
    plt.show()


# Define potentials
V0 = lambda x: 0
V1 = lambda x: k*x

def V2(x):
    if -3*h < x and x < 3*h:
        return k
    else:
        return 0

def V3(x):
    if x < -3*h:
        return -k
    elif -3*h < x and x < 3*h:
        return k*(-1 + 2 * ((x + 3*h)/(6*h)))
    else:
        return k
'''
beta_k_list = [0.1, 3, 5, 7, 10]
for elem in beta_k_list:
    beta_k = elem
    k = beta_k/beta
    print(k)
    plot_distribution(particles, V0)
    particles = np.zeros(N)
    plot_distribution(particles, V1)
    particles = np.zeros(N)
    plot_distribution(particles, V2)
    particles = np.zeros(N)
    plot_distribution(particles, V3)
    particles = np.zeros(N)
'''

''' Simple action potential below'''

L = 50
h = 1
V_Na0 = 1
V_K_0 = 1
Cc = 0.07*elemc #M/V
x_membrane = np.arange(-h, h)
C_p = 0.1 #mM - concentration per particle

N_Na = 60 + 1450 # inside + outside
N_K = 1400 + 50 # inside + outside
Na_pos = np.array([L/4]*1450 + [-L/4]*60)
K_pos = np.array([L/4]*50 + [-L/4]*1400)

steps = 400

# channel potentials

def V_channel(x, V0):
    if -h <= x and x <= h
        return V0
    else:
        return 0

#TODO: count # particles (Na and K) inside and outside of the membrane

def V_elec(x):
    # calculate Qc somewhere from distribution
    # (number of particles inside * C_p - outside concentration)
    return Qc/Cc














