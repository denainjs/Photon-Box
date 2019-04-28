## TODO: header explaining context
import qutip as q
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as random


## Simulation parameters
num_steps = 400 # number of steps in the simulation 
monte_carlo_steps = 10 # number of simulations in Monte-Carlo method

## Hilbert space parameters 
n_real = 20 # dimension of the Hilbert space describing the cavity 
n_estimate = 15 # maximum number of photon unambiguously measured
a_real = q.destroy(n_real) # annihilation operator; the creation operator is given by a.dag()
a_estimate = q.destroy(n_estimate) # annihilation operator; the creation operator is given by a.dag()

## Target
target = 3 # we want to stabilize the Fock state |target>
rho_target_real = q.fock_dm(n_real,target) # target density matrix
rho_target_estimate = q.fock_dm(n_estimate,target) # target density matrix

## Control parameters
c_1 = 1/(4*target + 2)
c_2 = 0.1
epsilon = 0.001 

## Cavity parameters
phi_bar_estimate =  np.pi/n_estimate # parameter for phi(n)
phi_bar_real =  np.pi/n_real # parameter for phi(n)
phi_R_estimate = np.pi/2 - phi_bar_estimate*target # parameter for phi(n)
phi_R_real = np.pi/2 - phi_bar_real*target # parameter for phi(n)

## Initial state
rho_init_real = q.coherent_dm(n_real,target**0.5) # initial density matrix
rho_init_estimate = q.coherent_dm(n_estimate,target**0.5) # initial density matrix


## Measurement operators
Mg_real = q.Qobj([[0 for i in range(n_real)] for j in range(n_real)])
Me_real = q.Qobj([[0 for i in range(n_real)] for j in range(n_real)])
Mg_estimate = q.Qobj([[0 for i in range(n_estimate)] for j in range(n_estimate)])
Me_estimate = q.Qobj([[0 for i in range(n_estimate)] for j in range(n_estimate)])
for i in range(0,n_real):
    Mg_real += np.cos((phi_R_real + phi_bar_real * i)/2) * q.fock_dm(n_real,i)
    Me_real += np.sin((phi_R_real + phi_bar_real * i)/2) * q.fock_dm(n_real,i)
for i in range(0,n_estimate):
    Mg_estimate += np.cos((phi_R_estimate + phi_bar_estimate * i)/2) * q.fock_dm(n_estimate,i)
    Me_estimate += np.sin((phi_R_estimate + phi_bar_estimate * i)/2) * q.fock_dm(n_estimate,i)

## Omission parameters
eta = 0.8 # probability of detecting the atom

def markov_chain(target):
    lyapunov_values_estimate = []
    lyapunov_values_real = []
    alpha_k_values = []
    rho_real = rho_init_real # actual (pure) state of the system
    rho_estimate = rho_init_estimate # estimated density over possible states of the system
    lyapunov_values_estimate.append(lyapunov(rho_estimate,rho_target_estimate))
    lyapunov_values_real.append(lyapunov(rho_real,rho_target_real))
    for i in range(num_steps):
        print("Iteration number:", i)
        # Measurement: project rho_real, project rho_estimate if measurement occurs
        proba_excited = np.real((Me_real*rho_real*Me_real.dag()).tr())
        mu = np.random.binomial(1, proba_excited) # 0 is ground, 1 is excited
        is_measured = np.random.binomial(1, eta)
        if mu == 0:
            rho_real = Mg_real*rho_real*Mg_real.dag()
            rho_real /= rho_real.tr()
            if is_measured == 1:
                rho_estimate = Mg_estimate*rho_estimate*Mg_estimate.dag()
                rho_estimate /= rho_estimate.tr()
        if mu == 1:
            rho_real = Me_real * rho_real * Me_real.dag()
            rho_real /= rho_real.tr()
            if is_measured == 1:
                rho_estimate = Me_estimate*rho_estimate*Me_estimate.dag()
                rho_estimate /= rho_estimate.tr()
        if is_measured == 0:
            rho_estimate = Mg_estimate*rho_estimate*Mg_estimate.dag() +                                                     Me_estimate*rho_estimate*Me_estimate.dag()
        # Control
        if np.real((rho_estimate * rho_target_estimate).tr()) >= epsilon:
            alpha_k = c_1 * (commutator(rho_target_estimate, a_estimate.dag() - a_estimate) * rho_estimate).tr()
        else:
            alpha_k = c_2 * np.sign(target - q.expect(a_estimate.dag() * a_estimate, rho_estimate))
        alpha_k_values.append(alpha_k)
        rho_real = q.displace(n_real, alpha_k) * rho_real * q.displace(n_real, -alpha_k)
        rho_estimate = q.displace(n_estimate, alpha_k) * rho_estimate * q.displace(n_estimate, -alpha_k)
        # Add the Lyapunov value to the list
        lyapunov_values_estimate.append(lyapunov(rho_estimate,rho_target_estimate))
        lyapunov_values_real.append(lyapunov(rho_real,rho_target_real)) 
    return lyapunov_values_real, lyapunov_values_estimate ,alpha_k_values

def commutator(op1, op2):
    return op1 * op2 - op2 * op1

def lyapunov(rho_1, rho_2): # in this case the Lyapunov function is the trace distance 
    return 1 - (rho_1 *rho_2).tr()

def plot_fidelity(lyapunov_values):
    steps = q.arange(num_steps + 1)
    plt.ylim(0.,1.)
    plt.plot(steps, lyapunov_values)
    plt.show()



def plot_many_lyapunov_values(num_curves):
    steps = q.arange(num_steps + 1)
    plt.ylim(0.,1.)
    real_average = np.zeros(num_steps + 1, dtype = 'complex128')
    estimate_average = np.zeros(num_steps + 1, dtype = 'complex128')
    for i in range(num_curves):
        markov_chain_results = markov_chain(target)
        lyapunov_values_real = markov_chain_results[0]
        lyapunov_values_estimate = markov_chain_results[1]
        real_average += lyapunov_values_real
        estimate_average += lyapunov_values_estimate
        plt.plot(steps, lyapunov_values_real, color = 'b')
        plt.plot(steps, lyapunov_values_estimate, color = 'r')
    plt.plot(steps, real_average/num_curves, color = 'g')
    plt.plot(steps, estimate_average/num_curves, color = 'orange')
    plt.show()

plot_many_lyapunov_values(5)

def plot_many_alpha_k(num_curves):
    steps = q.arange(num_steps)
    for i in range(num_curves):
        alpha_k_values = markov_chain(target)[1]
        plt.plot(steps, alpha_k_values)
    plt.show()

# plot_many_alpha_k(10)