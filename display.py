import matplotlib.pyplot as plt
import numpy as np

num_steps = 300

## Open loop

def display_open_loop(res):
    plt.xlabel("Iterations")
    plt.ylabel("Fidelity")
    plt.plot(np.arange(num_steps + 1), res)

## With control, without Lindblad, with the fidelity

def display_results_fidelity(res):
    fig2, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9,6))
    plt.rc('font',size=10)
    steps_plus = np.arange(num_steps + 1)
    steps = np.arange(num_steps)
    # Control
    control = res[2]
    ax1.set_xlabel("Iterations")
    ax1.set_ylabel("Control alpha")
    ax1.set_xlim(0, num_steps + 1)
    ax1.set_ylim(-0.1,0.1)
    ax1.plot(steps, control, color = 'black')
    # Real probabilities
    prob = res[0]
    ax2.set_xlabel("Iterations")
    ax2.set_ylabel("Real fidelity")
    ax2.set_xlim(0, num_steps + 1)
    ax2.set_ylim(0.0,1.0)
    ax2.plot(steps_plus, prob, color = 'b')
    # dist_estimate
    dist_estimate = 1 - np.array(res[1])
    ax3.set_xlabel("Iterations")
    ax3.set_ylabel("Estimated distance")
    ax3.set_xlim(0, num_steps + 1)
    ax3.set_ylim(0.0,1.0)
    ax3.plot(steps_plus, dist_estimate, color = 'c')
    plt.show()
    
def display_results_fidelity_lindblad(res):
    fig2, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9,6))
    plt.rc('font',size=10)
    steps_plus = np.arange(num_steps + 1)
    steps = np.arange(num_steps)
    # Control
    control = res[2]
    ax1.set_xlabel("Iterations")
    ax1.set_ylabel("Control alpha")
    ax1.set_xlim(0, num_steps + 1)
    ax1.set_ylim(-0.1,0.1)
    ax1.plot(steps, control, color = 'black')
    # Real probabilities
    prob = res[0]
    P_inf = res[3]
    P_sup = res[4]
    ax2.set_xlabel("Iterations")
    ax2.set_ylabel("Real fidelity")
    ax2.set_xlim(0, num_steps + 1)
    ax2.set_ylim(0.0,1.0)
    ax2.plot(steps_plus, prob, color = 'b')
    ax2.plot(steps_plus, P_inf, color = 'r')
    ax2.plot(steps_plus, P_sup, color = 'g')
    # dist_estimate
    dist_estimate = res[1]
    ax3.set_xlabel("Iterations")
    ax3.set_ylabel("Distance 1-W0")
    ax3.set_xlim(0, num_steps + 1)
    ax3.set_ylim(0.0,1.0)
    ax3.plot(steps_plus, dist_estimate, color = 'c')
    plt.show()

# and display_results_W
