## TODO: header explaining context
import qutip as q
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as random

num_steps = 200 # number of steps in the simulation 
monte_carlo_steps = 50 # number of simulations in Monte-Carlo method

N = 10 # dimension of the Hilbert space describing the cavity 
n_max = 10 # maximum number of photon unambiguously measured

target = 3 # we want to stabilize the Fock state |target>
alpha = target**0.5 # the initial coherent state

phi_bar =  np.pi/n_max # parameter for phi(n)
phi_R = np.pi/2 - phi_bar*target # parameter for phi(n)

a = q.destroy(N) # annihilation operator; the creation operator is given by a.dag()

rho_target = q.fock_dm(N,target) # target density matrix
# print("Target density matrix: \n", rho_target)

rho_init = q.coherent_dm(N,alpha) # initial density matrix
# print("Initial density matrix: \n", rho_init)
# print("Its trace is:", rho_init.tr())
# print("The initial mean number of photons is:", q.expect(a.dag() * a ,rho_init))

## TODO: define the measurement operators
Mg = q.Qobj([[0 for i in range(N)] for j in range(N)])
Me = q.Qobj([[0 for i in range(N)] for j in range(N)])
for i in range(0,N):
    Mg += np.cos((phi_R + phi_bar * i)/2) * q.fock_dm(N,i)
    Me += np.sin((phi_R + phi_bar * i)/2) * q.fock_dm(N,i)
# print("Measurement of the ground state: \n", Mg)
# print("Measurement of the excited state: \n", Me)


def markov_chain(target):
    lyapunov_values = []
    rho = rho_init
    lyapunov_values.append(lyapunov(rho,rho_target))
    for i in range(num_steps):
        proba_excited = np.real((Me*rho*Me.dag()).tr())
        y = np.random.binomial(1, proba_excited) # 0 is ground, 1 is excited
        if y == 0:
            # print("Atom measured in Ground state: \n")
            rho = Mg*rho*Mg.dag()
            rho /= rho.tr()
        if y == 1:
            # print("Atom measured in Excited state: \n")
            rho = Me*rho*Me.dag()
            rho /= rho.tr()
        lyapunov_values.append(lyapunov(rho,rho_target)) 
        # print("Iteration ", i , "the density matrix is \n", rho)
        # print("Its trace is:", rho.tr())
        # print("The mean value is:", q.expect(a.dag() * a, rho))
    return lyapunov_values


def lyapunov(rho_1, rho_2): # in this case the Lyapunov function is the trace distance 
    return 1 - (rho_1 *rho_2).tr()

def plot_fidelity(lyapunov_values):
    steps = q.arange(num_steps + 1)
    plt.ylim(0.,1.)
    plt.plot(steps, lyapunov_values)
    plt.show()

def monte_carlo_average_fidelity(target, monte_carlo_steps, num_steps):
    sum_lyapunov_values = np.array([0 for i in range(num_steps+1)], dtype = "complex128")
    for _ in range(monte_carlo_steps):
        sum_lyapunov_values += np.array(markov_chain(target))
    return sum_lyapunov_values/monte_carlo_steps
        


def plot_many_lyapunov_values(num_curves):
    steps = q.arange(num_steps + 1)
    plt.ylim(0.,1.)
    average = np.array([0 for i in range(num_steps+1)],dtype = "complex128")
    for i in range(num_curves):
        lyapunov_values = markov_chain(target)
        average += np.array(lyapunov_values)
        plt.plot(steps, lyapunov_values)
    average /= num_curves
    plt.plot(steps, average)
    plt.show()

plot_many_lyapunov_values(50)
