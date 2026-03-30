import numpy as np

# Bustamante's objective function.
# According to the parallelization results in output_dict performed at a particular force,
# it calculates the objective (squared error) function between expected number of cycles
# calculated from the simulations and the kinetic model derived from Bustamante's single molecule experiment
#
def Bustamante_objective(output_dict, force, circuit_length=2200):

    # Initial conditions
    Lk0 = circuit_length / 10.5 # Relaxed linking number
    total_rotations_induced = [] # This is a list that contains a distribution of induced rotaitons per binding event

    for i, system_output in enumerate(output_dict):

        n_sims = len(system_output['sites_df_list'])

        for n_sim in range(n_sims):
            sites_df = system_output['sites_df_list'][n_sim]
            # enzymes_df = system_output['enzymes_df_list'][n_sim]
            total_rotations_induced += rotations_induced_per_binding_event(sites_df, Lk0)

    # Calculate objective
    n_cycles_sim = -np.mean(np.array(total_rotations_induced) / 2.0) # Expected number of cycles from simulation data
    p_cycles_exp, n_cycles_exp = mechano_gyrase(force)  # Expected number of cycles from single molecule experiments

    error = (n_cycles_sim - n_cycles_exp)**2 # Computes the error between

    return error, total_rotations_induced

# Given a sites_df, it computes the rotations induced per binding event. This function is intended to be
# applied to topoisomerases.
def rotations_induced_per_binding_event(sites_df, Lk0):
    mask = sites_df['type'] == 'circuit'
    nenzymes = sites_df[mask]['#enzymes'].to_numpy()  #
    sigma = sites_df[mask]['superhelical'].to_numpy()
    dLk = sigma * Lk0  # This is the linking difference
    dif = np.diff(dLk[1:])  # Differences

    diff = np.diff(nenzymes[1:])
    p_start = np.where(diff == 1)[0] + 1
    p_end = np.where(diff == -1)[0] + 1

    rotations_induced = []
    for start, end in zip(p_start, p_end):
        twist = dif[start:end].sum()  # Summing the duration
        rotations_induced.append(twist)
    return rotations_induced

# Kinetic model derived in Bustamante's experiment. It calculates the probability of performing an additional cycle
# and the expected number of cycles at a force F.
def mechano_gyrase(F, kwrap_koff=85.0, delta_xt=31.0):
    B = F * delta_xt /4.1 # delta x = 31, KBT = 4.1
    C = kwrap_koff*np.exp(-B)
    p_cycles = C/(C+1)
    n_cycles = 1/(1-p_cycles)
    return p_cycles, n_cycles # Average number of cycles

