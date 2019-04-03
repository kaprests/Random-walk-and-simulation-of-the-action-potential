#TODO:
# Fix error(s)
# Proposal: make some simple test potentials (time independent), and make sure everything common with earlier
# tasks is functioning.


import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
from scipy.constants import elementary_charge as elemc, Boltzmann as kb

h = 1 # step length
L = 50 # h. length of system
T = 273 + 37 # K, temperature
betaV0 = 0 # beta * V0, beta = 1/kb*T, V0 some constant and kb = Boltzmanns constant
beta = 1/(kb*T)
beta_k = 1
Cc = 0.07*elemc*1e3 #CmM/V
Qc_out = 150 #mM
C_p = 0.1 # mM
N_Na = 50+1450
N_K = 1400+50
Na_pos = np.array([-L//4]*50 + [L//4]*1450)
K_pos = np.array([-L//4]*1400 + [L//4]*50)
steps = 500
pos_vec = np.arange(-steps, steps+1)


# potentials

# time independent potential of the channels
def V_channel(x, V_0):
    if -h <= x and x <= h:
        return V_0
    else:
        return 0


V_channel = np.vectorize(V_channel)


# Time dependent potential
def V_elec(Na_pos, K_pos):
    Na_in = Na_pos[Na_pos < -h].size
    K_in = K_pos[K_pos < -h].size
    Qc_in = (Na_in + K_in)*C_p
    Qc = Qc_in - Qc_out
    return elemc*Qc/Cc # volts


# Returns probability of a single particle stepping to the left
def P_min(x, V_vec):
    V1 = V_vec[steps + int(x) -h]
    V2 = V_vec[steps + int(x) +h]
    rel_prob = np.exp(-beta_k*(V1 - V2))
    Pp = 1/(1 + rel_prob)
    Pm = 1 - Pp
    if V1 - V2 != 0: # debug prints :'((
        print("")
        print(x)
        print(V1 - V2)
        print(V1)
        print(V2)
        print(V_vec[steps-3: steps+4])
        print(Pm)
    return Pm


# Function looping through time and performing the simulation. Returns a vector with the
# time independent potential values at each time step.
def rand_walk(Na_pos_vec, K_pos_vec, V_Na_vec, V_K_vec, p_vec, V_el_func, P_min_func):
    # array to store values of the time dependent potential
    Ve_vec = np.zeros(steps)

    # Loop over each time step
    for i in range(steps):
        # Calculate time dependent potential for the current timestep and store it in Ve_Vec
        Ve = V_el_func(Na_pos_vec, K_pos_vec)
        Ve_vec[i] = Ve

        # add Ve to the time independent potential vectors to get total potential
        V_Na_tot = V_Na_vec + np.heaviside(p_vec, 0.5)*Ve*elemc #Ve*elemc shoud be joules, but who knows :)
        V_K_tot = V_K_vec + np.heaviside(p_vec, 0.5)*Ve*elemc
    
        #print(V_Na_tot[steps-3:steps+3])

        # Generate vectors with 1 and -1 (vector containing the next step for each particle)
        steps_vec_Na = np.random.rand(Na_pos_vec.size) # vector with probabilities
        steps_vec_K = np.random.rand(K_pos_vec.size) # vector with probabilities
        # loop over all probabilities and determine if it should be a step to the right or left
        for j in range(steps_vec_Na.size):
            if steps_vec_Na[j] >= P_min_func(Na_pos_vec[j], V_Na_tot):
                steps_vec_Na[j] = 1
        for j in range(steps_vec_K.size):
            if steps_vec_K[j] >= P_min_func(K_pos_vec[j], V_K_tot):
                steps_vec_K[j] = 1
        # makes every step thats not a right step become a left step, maybe cleaner to use else statements for this
        steps_vec_Na[steps_vec_Na != 1] = -1
        steps_vec_K[steps_vec_K != 1] = -1

        # Have every particle take one step according to the step vector
        Na_pos_vec += steps_vec_Na.astype(int)
        K_pos_vec += steps_vec_K.astype(int)

        # boundaries at +- L/2, if a particle is on the boundary it has to move one step into the system
        Na_pos_vec[Na_pos_vec < -L/2] = -L/2 + 1
        Na_pos_vec[Na_pos_vec > L/2] = L/2 - 1
        K_pos_vec[K_pos_vec < -L/2] = -L/2 + 1
        K_pos_vec[K_pos_vec > L/2] = L/2 - 1
    return Ve_vec


# Plotting function
def plot_dist(Na_pos_vec, K_pos_vec, V_Na_vec, V_k_vec, p_vec, V_el_func, P_min_func, t_vec, V_name="hm"):
    Ve_vec = rand_walk(Na_pos_vec, K_pos_vec, V_Na_vec, V_k_vec, p_vec,  V_el_func, P_min_func)

    plt.hist(Na_pos_vec, bins="auto")
    plt.show()
    plt.hist(K_pos_vec, bins="auto")
    plt.show()
    plt.plot(t_vec, Ve_vec)
    plt.show()


V_Na = V_channel(pos_vec, betaV0/beta)
V_K = V_channel(pos_vec, betaV0/beta)
timesteps = np.arange(steps)

print(V_Na)
print(V_K)

plot_dist(Na_pos, K_pos, V_Na, V_K, pos_vec, V_elec, P_min, timesteps)



