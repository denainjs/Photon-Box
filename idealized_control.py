## TODO: header explaining context

import numpy as np 
import numpy.random as random
import matplotlib.pyplot as plt

num_steps = 50 # number of steps in the simulation 
monte_carlo_steps = 100 # number of simulations in Monte-Carlo method

N = 10 # dimension of the Hilbert space describing the cavity 
n_max = 10 # maximum number of photon unambiguously measured

target = 3 # we want to stabilize the Fock state |target>
alpha = np.sqrt(target) # we start with the coherent state with mean number of photons |alpha|^2 

phi_bar =  np.pi/n_max # parameter for phi(n)
phi_R = np.pi/2 - phi_bar*target # parameter for phi(n)

rho_target = np.zeros((N,N))
rho_target[target][target] = 1

Mg = np.zeros((N,N))
Me = np.zeros((N,N))



J = np.zeros((N,N))
for i in range(N):
    J[i][i] = i

for i in range(N):
    Mg[i][i] = np.cos(phi_R + phi_bar * i)
    Me[i][i] = np.sin(phi_R + phi_bar * i)
    
# print("The matrix Mg is: \n", Mg)
# print()
# print("The matrix Me is: \n", Me)
# print()

def initialize_rho(alpha): # creates the density matrix corresponding to the pure coherent state |alpha>
    rho = np.zeros((N,N))
    value = np.exp(-np.abs(alpha)**2)
    rho[0][0] = value
    for i in range(1,N):
        value *=  np.abs(alpha)**2/i
        rho[i][i] = value
    return rho

## TODO: from transpose to adjoint (requires to handle numpy.matrix)
# for now, just tranpose



def markov_chain(alpha):
    fidelities = []
    rho = initialize_rho(alpha) 
    fidelities.append(fidelity(rho,rho_target))
    # print("Initial density matrix: \n", rho)
    # print("its trace is:", np.trace(rho))
    # print()
    # print("the mean value is:", np.trace(np.dot(J,rho)))
    # print()
    for i in range(num_steps):
        y = np.random.binomial(1, np.trace(np.dot(np.dot(Me,rho),np.transpose(Me)))) # 0 is ground, 1 is excited
        if y == 0:
            # print("ground: \n")
            rho = np.dot(np.dot(Mg,rho),np.transpose(Mg))
            rho /= np.trace(rho)
        if y == 1:
            # print("excited: \n")
            rho = np.dot(np.dot(Me,rho),np.transpose(Me))
            rho /= np.trace(rho)
        fidelities.append(fidelity(rho,rho_target)) 
        # print("iteration ", i , "the density matrix is \n")
        # print()
        # print("the density matrix is: \n", rho)
        # print()
        # print("its trace is:", np.trace(rho))
        # print()
        # print("the mean value is:", np.trace(np.dot(J,rho)))
        # print()
    return fidelities


def fidelity(rho_1, rho_2):
    return 1 - np.trace(np.dot(rho_1, rho_2))

def plot_fidelity(fidelities):
    steps = np.arange(num_steps + 1)
    plt.ylim(0.,1.)
    plt.plot(steps, fidelities)
    plt.show()

def monte_carlo_average_fidelity(alpha, monte_carlo_steps, num_steps):
    sum_fidelities = np.array([0 for i in range(num_steps+1)],dtype = "float64")
    for _ in range(monte_carlo_steps):
        sum_fidelities += np.array(markov_chain(alpha))
    return sum_fidelities/monte_carlo_steps
        



plot_fidelity(monte_carlo_average_fidelity(alpha, monte_carlo_steps, num_steps))