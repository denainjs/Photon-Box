import numpy as np
import numpy.random as random
import qutip as q

## All the parameters

num_steps = 300

n_real = 20 # dimension of the Hilbert space describing the cavity 
n_estimate = 15 # maximum number of photon unambiguously measured

target = 3 # we want to stabilize the Fock state |target>
rho_target_real = q.fock_dm(n_real,target) # target density matrix
rho_target_estimate = q.fock_dm(n_estimate,target) # target density matrix

phi_bar_estimate =  np.pi/n_estimate # parameter for phi(n)
phi_bar_real =  np.pi/n_real # parameter for phi(n)
phi_R_estimate = np.pi/2 - phi_bar_estimate*target # parameter for phi(n)
phi_R_real = np.pi/2 - phi_bar_real*target # parameter for phi(n)

rho_init_real = q.coherent_dm(n_real,target**0.5) # initial density matrix
rho_init_estimate = q.coherent_dm(n_estimate,target**0.5) # initial density matrix

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
    
c_1 = 1/(4*target + 2)
c_2 = 0.1
epsilon = 0.001 
a_real = q.destroy(n_real) # annihilation operator; the creation operator is given by a_real.dag()
a_estimate = q.destroy(n_estimate) # annihilation operator; the creation operator is given by a_estimate.dag()


eta = 0.8 

## Useful functions

def displacement(rho,alpha, dimension):
    return q.displace(dimension, alpha)*rho*q.displace(dimension, -alpha)

def commutator(op1, op2):
    return op1 * op2 - op2 * op1

def fidelity(rho, rho_target):
    return (rho * rho_target).tr()

## Kraus map functions

def measure_projection(rho, measurement_operator):
    res = measurement_operator * rho * measurement_operator.dag()
    return res/res.tr()

def kraus_map_omission(rho, Mg_estimate, Me_estimate):
    return Mg_estimate*rho*Mg_estimate.dag() + Me_estimate*rho*Me_estimate.dag()

## The actual Markov Chain

def omission_markov_chain():
    fidelity_values_estimate = []
    fidelity_values_real = []
    alpha_k_values = []
    rho_real = rho_init_real # actual (pure) state of the system
    rho_estimate = rho_init_estimate # estimated density over possible states of the system
    fidelity_values_estimate.append(fidelity(rho_estimate,rho_target_estimate))
    fidelity_values_real.append(fidelity(rho_real,rho_target_real))
    for i in range(num_steps):
        # Measurement: project rho_real and rho_estimate
        proba_excited = np.real((Me_real*rho_real*Me_real.dag()).tr())
        mu = np.random.binomial(1, proba_excited) # 0 is ground, 1 is excited
        is_measured = np.random.binomial(1, eta)
        if mu == 0:
            rho_real = measure_projection(rho_real, Mg_real)
            if is_measured == 1:
                rho_estimate = measure_projection(rho_estimate, Mg_estimate)
        if mu == 1:
            rho_real = measure_projection(rho_real, Me_real)
            if is_measured == 1:
                rho_estimate = measure_projection(rho_estimate, Me_estimate)
        if is_measured == 0:
            rho_estimate = kraus_map_omission(rho_estimate, Mg_estimate, Me_estimate)
        # Control
        if np.real((rho_estimate * rho_target_estimate).tr()) >= epsilon:
            alpha_k = c_1 * (commutator(rho_target_estimate, a_estimate.dag() - a_estimate) * rho_estimate).tr()
        else:
            alpha_k = c_2 * np.sign(target - q.expect(a_estimate.dag() * a_estimate, rho_estimate))
        # Add the value of the control to the list
        alpha_k_values.append(alpha_k)
        # Apply the control
        rho_real = displacement(rho_real,alpha_k, n_real)
        rho_estimate = displacement(rho_estimate,alpha_k, n_estimate)
        # Add the value of the real and estimated fidelities to the list
        fidelity_values_real.append(fidelity(rho_real,rho_target_real)) 
        fidelity_values_estimate.append(fidelity(rho_estimate,rho_target_estimate))
    return fidelity_values_real, fidelity_values_estimate, alpha_k_values