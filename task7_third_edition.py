import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import Boltzmann as kb, elementary_charge as elemc


h = 1 # step length
L = 50 # system length
T = 273 + 37
beta_V0_Na = 3
beta_V0_K = 3
beta = 1/(kb*T)
Cc = 0.07*elemc*1e3 # mMC/V
Qc_out = 150 # mM
Cp = 0.1 # consentratioin per particle
N_Na = 50 + 1450
N_K = 1400+50
Na_pos = np.array([-L//2]*50 + [L//2]*1450)
K_pos = np.array([-L//2]*1400 + [L//2]*50)
t_steps = 500
position_vec = np.arange(-t_steps, t_steps +1)


def V_channel(x, V_0):
    if -h <= x and x <= h:
        return V_0
    else:
        return 0
V_channel = np.vectorize(V_channel)


def V_electric(Na_pos_vec, K_pos_vec):
    num_Na_in = Na_pos_vec[Na_pos_vec < -h].size
    num_K_in = K_pos_vec[K_pos_vec < -h].size
    Qc_in = (num_Na_in + num_K_in)*Cp
    Qc = Qc_in - Qc_out
    return elemc*Qc/Cc # Volts


def P_min(x, V_vec):
    V1 = V_vec[t_steps + int(x) -h]
    V2 = V_vec[t_steps + int(x) +h]
    rel_prob = np.exp(-beta*(V1 - V2))
    Pp = 1/(1 + rel_prob)
    Pm = 1- Pp
    return Pm


def rand_walk(Na_pos_vec, K_pos_vec, V_Na_vec, V_K_vec, pos_vec, P_min_func, V_elec_func):
    Ve_vec = np.zeros(t_steps)
    for i in range(t_steps):
        Ve = V_elec_func(Na_pos_vec, K_pos_vec)
        Ve_vec[i] = Ve

        V_Na_tot = V_Na_vec + np.heaviside(-pos_vec, 0.5)*Ve*elemc
        V_K_tot = V_K_vec + np.heaviside(-pos_vec, 0.5)*Ve*elemc

        Na_steps_vec = np.random.rand(Na_pos_vec.size)
        K_steps_vec = np.random.rand(K_pos_vec.size)

        for i in range(Na_steps_vec.size):
            if Na_steps_vec[i] >= P_min_func(Na_pos_vec[i], V_Na_tot):
                Na_steps_vec[i] = 1
            else:
                Na_steps_vec[i] = -1

        for i in range(K_steps_vec.size):
            if K_steps_vec[i] >= P_min_func(K_pos_vec[i], V_Na_tot):
                K_steps_vec[i] = 1
            else:
                K_steps_vec[i] = -1

        Na_pos_vec += Na_steps_vec.astype(int)
        K_pos_vec += K_steps_vec.astype(int)

        Na_pos_vec[Na_pos_vec < -L/2] = -L/2 +1
        Na_pos_vec[Na_pos_vec > L/2] = L/2 -1

        K_pos_vec[K_pos_vec < -L/2] = -L/2 +1
        K_pos_vec[K_pos_vec > L/2] = L/2 -1

    return Ve_vec


def plot_distribution(Na_pos_vec, K_pos_vec, V_Na_vec, V_K_vec, pos_vec, t_vec, P_min_func, V_elec_func):
    Ve_vec = rand_walk(Na_pos_vec, K_pos_vec, V_Na_vec, V_K_vec, pos_vec, P_min_func, V_elec_func)

    plt.hist(Na_pos, bins="auto")
    plt.show()
    plt.hist(K_pos, bins="auto")
    plt.show()
    plt.plot(t_vec, Ve_vec)
    plt.show()


V_Na = V_channel(position_vec, beta_V0_Na/beta)
V_K = V_channel(position_vec, beta_V0_K/beta)
time_vec = np.arange(t_steps)


plot_distribution(Na_pos, K_pos, V_Na, V_K, position_vec, time_vec, P_min, V_electric)


