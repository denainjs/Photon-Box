## TODO: header explaining context
import qutip as q
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as random

## Simulation parameters
num_steps = 200 # number of steps in the simulation 
monte_carlo_steps = 50 # number of simulations in Monte-Carlo method

## Hilbert space parameters 
N = 10 # dimension of the Hilbert space describing the cavity 
n_max = 10 # maximum number of photon unambiguously measured
a = q.destroy(N) # annihilation operator; the creation operator is given by a.dag()

## Target
target = 3 # we want to stabilize the Fock state |target>
rho_target = q.fock_dm(N,target) # target density matrix
# print("Target density matrix: \n", rho_target)

## Control parameters
c_1 = 1/(4*target + 2)
c_2 = 0.1
epsilon = 0.001 

## Cavity parameters
phi_bar =  np.pi/n_max # parameter for phi(n)
phi_R = np.pi/2 - phi_bar*target # parameter for phi(n)

## Initial state
rho_init = q.coherent_dm(N,target**0.5) # initial density matrix
# print("Initial density matrix: \n", rho_init)
# print("Its trace is:", rho_init.tr())
# print("The initial mean number of photons is:", q.expect(a.dag() * a ,rho_init))


## Measurement operators
Mg = q.Qobj([[0 for i in range(N)] for j in range(N)])
Me = q.Qobj([[0 for i in range(N)] for j in range(N)])
for i in range(0,N):
    Mg += np.cos((phi_R + phi_bar * i)/2) * q.fock_dm(N,i)
    Me += np.sin((phi_R + phi_bar * i)/2) * q.fock_dm(N,i)
# print("Measurement of the ground state: \n", Mg)
# print("Measurement of the excited state: \n", Me)

def markov_chain(target):
    lyapunov_values = []
    alpha_k_values = []
    rho = rho_init
    lyapunov_values.append(lyapunov(rho,rho_target))
    for i in range(num_steps):
        # Measurement
        proba_excited = np.real((Me*rho*Me.dag()).tr())
        y = np.random.binomial(1, proba_excited) # 0 is ground, 1 is excited
        if y == 0:
            # print("Atom measured in Ground state: \n")
            rho = Mg*rho*Mg.dag()
            rho /= rho.tr()
        if y == 1:
            # print("Atom measured in Excited state: \n")
            rho = Me * rho * Me.dag()
            rho /= rho.tr()
        # Control
        if np.real((rho * rho_target).tr()) >= epsilon:
            alpha_k = c_1 * (commutator(rho_target, a.dag() - a) * rho).tr()
        else:
            alpha_k = c_2 * np.sign(target - q.expect(a.dag() * a, rho))
        alpha_k_values.append(alpha_k)
        rho = q.displace(N, alpha_k) * rho * q.displace(N, -alpha_k)
        # Add the Lyapunov value to the list
        lyapunov_values.append(lyapunov(rho,rho_target)) 
        # print("Iteration ", i , "the density matrix is \n", rho)
        # print("Its trace is:", rho.tr())
        # print("The mean value is:", q.expect(a.dag() * a, rho))
    return lyapunov_values,alpha_k_values

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
    for i in range(num_curves):
        lyapunov_values = markov_chain(target)[0]
        plt.plot(steps, lyapunov_values)
    plt.show()

# plot_many_lyapunov_values(10)

def plot_many_alpha_k(num_curves):
    steps = q.arange(num_steps)
    for i in range(num_curves):
        alpha_k_values = markov_chain(target)[1]
        plt.plot(steps, alpha_k_values)
    plt.show()

plot_many_alpha_k(10)