import numpy as np 
import numpy.random as random

num_steps = 100 # number of steps in the simulation 
MC = 100 # number of simulations in Monte-Carlo method

N = 10 # dimension of the Hilbert space describing the cavity 
n_max = 10 # maximum number of photon unambiguously measured

target = 3 # we want to stabilize the Fock state |target>
alpha = np.sqrt(target) # we start with the coherent state with mean number of photons |alpha|^2 

phi_bar =  np.pi/n_max # parameter for phi(n)
phi_R = np.pi/2 - phi_bar*target # parameter for phi(n)

J = np.zeros((N,N)) # matrix of the operator N in the Fock basis
for i in range(N):
    J[i][i] = i

Mg = np.cos(phi_R + phi_bar * J) # quantum measurement operator, ground state
Me = np.sin(phi_R + phi_bar * J) # quantum measurement operators, excited state



def initialize_rho(alpha): # creates the density matrix corresponding to the pure coherent state |alpha>
    rho = np.zeros((N,N))
    value = np.exp(-np.abs(alpha)**2)
    rho[0][0] = value
    for i in range(1,N):
        value *=  np.abs(alpha)**2/i
        rho[i][i] = value
    return rho

## TODO: from transpose to adjoint (requires to handle numpy.matrix)

def markov_chain(alpha):
    rho = initialize_rho(alpha) 
    print("Initial density matrix:", rho)
    for i in range(num_steps):
        y = np.random.randint(2) # 0 is ground, 1 is excited
        if y == 0:
            rho = np.dot(np.dot(Mg,rho),np.transpose(Mg))
            rho /= np.trace(np.dot(np.dot(Mg,rho),np.transpose(Mg)))
        if y == 1:
            rho = np.dot(np.dot(Me,rho),np.transpose(Me))
            rho /= np.trace(np.dot(np.dot(Me,rho),np.transpose(Me)))
        print()
        print("iteration ", i , "the density matrix is \n")
        print()
        print(rho)
        print()

def fidelity(rho_1, rho_2):
    return 1 - np.trace(np.dot(rho_1, rho_2))

markov_chain(alpha)