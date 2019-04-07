import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
from scipy.constants import elementary_charge as elemc, Boltzmann as kb

h = 1 # step length
L = 50 # h. length of system
T = 273 + 37 # K, temperature
betaV0_Na = 0 # beta * V0_Na, beta = 1/kb*T, V0 some constant and kb = Boltzmanns constant
betaV0_K = 0 # beta * V0_K, beta = 1/kb*T, V0 some constant and kb = Boltzmanns constant
beta = 1/(kb*T)
Cc = 0.07*elemc*1e3 #CmMC/V
Qc_out = 150 #mM
C_p = 0.1 # mM
N_Na = 50 + 1450
N_K = 1400 + 50
Na_pos = np.array([-L//4]*50 + [L//4]*1450)
K_pos = np.array([-L//4]*1400 + [L//4]*50)
steps = 500
pos_vec = np.arange(-steps, steps+1)

# Regulates potential according to voltage difference
min_vol = -70 #mV
max_vol = 30 #mV
open_pot = 1


# potentials
def V_channel(pos_vec, V_0):
    channel_vec = np.array([V_0]*(2*h +1))
    zero_vec = np.zeros((pos_vec.size - channel_vec.size)//2)
    V_vec = np.concatenate((zero_vec, channel_vec))
    V_vec = np.concatenate((V_vec, zero_vec))
    return V_vec


# Time dependent potential
def V_elec(Na_pos, K_pos):
    Na_in = Na_pos[Na_pos < -h].size
    K_in = K_pos[K_pos < -h].size
    Qc_in = (Na_in + K_in)*C_p
    Qc = Qc_in - Qc_out
    return elemc*Qc/Cc #volts


# Returns probability of a single particle stepping to the left
def P_min(x, V_vec):
    V1 = V_vec[steps + int(x) -h]
    V2 = V_vec[steps + int(x) +h]
    rel_prob = np.exp(-beta*(V1 - V2))
    Pp = 1/(1 + rel_prob)
    Pm = 1 - Pp
    return Pm


def volt_reg(current_voltage):
    if current_voltage <= min_vol:
        Na_pot = V_channel(pos_vec, open_pot/beta)
        K_pot = V_channel(pos_vec, betaV0_K/beta)
    if current_voltage >= max_vol:
        Na_pot = V_channel(pos_vec, betaV0_Na/beta)
        K_pot = V_channel(pos_vec, open_pot/beta)
    return Na_pot, K_pot


# Pumps ions in and out membrane, step counter defined in rand_walk()
def pump(Na_pos, K_pos):
    N_Na_in = Na_pos[Na_pos < -h].size
    N_K_out = K_pos[K_pos > h].size
    if N_Na_in <= 3 or N_K_out <= 2:
        return None
    Na_pos.sort()
    K_pos.sort()
    Na_index = np.argwhere(Na_pos < -h)
    K_index = np.argwhere(K_pos > h)
    Na_pos[Na_index[-3:]] = h
    K_pos[K_index[:2]] = -h
    return Na_pos, K_pos


# Function looping through time and performing the simulation. Returns a vector with the
# time dependent potential values at each time step.
def rand_walk(Na_pos_vec, K_pos_vec, V_Na_vec, V_K_vec, p_vec, V_el_func, P_min_func):
    # array to store values of the time dependent potential
    Ve_vec = np.zeros(steps)

    # Loop over each time step
    for i in range(steps):
        # Calculate time dependent potential for the current timestep and store it in Ve_Vec
        Ve = V_el_func(Na_pos_vec, K_pos_vec)
        Ve_vec[i] = Ve

        # add Ve to the time independent potential vectors to get total potential
        V_Na_tot = V_Na_vec + np.heaviside(-p_vec, 0.5)*Ve*elemc
        V_K_tot = V_K_vec + np.heaviside(-p_vec, 0.5)*Ve*elemc

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


V_Na = V_channel(pos_vec, betaV0_Na/beta)
V_K = V_channel(pos_vec, betaV0_K/beta)
timesteps = np.arange(steps)

plot_dist(Na_pos, K_pos, V_Na, V_K, pos_vec, V_elec, P_min, timesteps)

