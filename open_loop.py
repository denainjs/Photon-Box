import numpy as np
import numpy.random as random
import qutip as q

## All the parameters

num_steps = 300

n = 10

target = 3 # we want to stabilize the Fock state |target>
rho_target = q.fock_dm(n,target)

rho_init = q.coherent_dm(n,target**0.5)

phi_bar =  np.pi/n 
phi_R = np.pi/2 - phi_bar*target 

Mg = q.Qobj([[0 for i in range(n)] for j in range(n)])
Me = q.Qobj([[0 for i in range(n)] for j in range(n)])
for i in range(0,n):
    Mg += np.cos((phi_R + phi_bar * i)/2) * q.fock_dm(n,i)
    Me += np.sin((phi_R + phi_bar * i)/2) * q.fock_dm(n,i)



## Useful functions

def fidelity(rho, rho_target):
    return (rho * rho_target).tr()

## Kraus map functions

def measure_projection(rho, measurement_operator):
    res = measurement_operator * rho * measurement_operator.dag()
    return res/res.tr()

## The actual Markov Chain

def open_loop_markov_chain(target):
    fidelity_values = []
    rho = rho_init
    fidelity_values.append(np.real(fidelity(rho,rho_target))) # fidelity is a real number
    for i in range(num_steps):
        proba_excited = np.real((Me*rho*Me.dag()).tr())
        y = np.random.binomial(1, proba_excited) # 0 is ground, 1 is excited
        if y == 0:
            rho = measure_projection(rho, Mg)
        if y == 1:
            rho = measure_projection(rho, Me)
        fidelity_values.append(np.real(fidelity(rho,rho_target))) # fidelity is a real number 
    return fidelity_values