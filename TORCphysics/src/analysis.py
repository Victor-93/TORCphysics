import numpy as np
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em


# This function inputs x, which indicates the number of enzymes bound to the corresponding site at time k
def build_1_signal(x):
    frames = len(x)
    signal = np.zeros(frames)
    for k in range(frames):
        if x[k] > 0:
            signal[k] = 1
    return signal


# Build signals by site type. Returns list of signals and list of site names (with same order than signals).
def build_signal_by_type(sites_df, my_type):
    mask = sites_df['type'] == my_type
    my_df = sites_df[mask]
    sites_names = my_df.drop_duplicates(subset='name')['name']
    signals = []
    names = []
    for name in sites_names:
        mask = my_df['name'] == name
        signal = my_df[mask]['#enzymes'].to_numpy()
        signals.append(build_1_signal(signal))
        names.append(name)
    return signals, names


# Build all signals in input sites_df. Returns list of signals and list of names
def build_signals(sites_df):
    sites_names = sites_df.drop_duplicates(subset='name')['name']
    signals = []
    names = []
    for name in sites_names:
        mask = sites_df['name'] == name
        signal = sites_df[mask]['#enzymes'].to_numpy()
        signals.append(build_1_signal(signal))
        names.append(name)
    return signals, names


# This one computes the cross-correlation hyper-matrix.
# Input: You give it a list of signals (each entry is one signal) that you want to calculate the cross-correlation, and
#        also give the timestep dt.
# Output: 1.- Returns a hyper-matrix (numpy) in the form [m=signal index,n=signal index,cross-correlations],
#             where diagonal elements are auto-correlations and off-diagonals are cross-correlations.
#         2.- It also returns the lag, which is the
def cross_correlation_hmatrix(signals, dt):
    if len(signals) <= 0:
        print('You gave me an empty list of signals, what is wrong with you?')
        return None, None
    frames = len(signals[0])
    n = len(signals)
    autocorr = np.zeros((n, frames))  # This is the hyper-matrix
    matrix = np.zeros((n, n, frames))  # This is the hyper-matrix
    lag = np.arange(-frames * dt * .5, frames * dt * .5, dt)

    # Let's first compute the auto-correlations
    for i, signal_i in enumerate(signals):
        autocorr[i, :] = np.correlate(signal_i, signal_i, "same")

    # Then the cross-correlation
    for i, signal_i in enumerate(signals):
        for j, signal_j in enumerate(signals):
            matrix[i, j, :] = np.correlate(signal_i, signal_j, 'same')
            if np.max(matrix[i, j, :]) > 0:
                matrix[i, j, :] = matrix[i, j, :] / np.sqrt(np.max(autocorr[i, :]) * np.max(autocorr[j, :]))

    return matrix, lag


# It computes the activity curves of site as a function of the supercoiling density.
# Input: 1.- Site. Optionals= sigma_min, sigma_max and delta_sigma
# Output: 1.- Modulated rates. 2.- Supercoiling density.
def site_activity_curves(site, sigma_min=-.2, sigma_max=.1, delta_sigma=.001):
    sigma = np.arange(sigma_min, sigma_max, delta_sigma)

    # Simple poisson process (constant binding)
    if site.site_model == 'poisson' or site.site_model == 'Poisson':
        rate = site.k_min * np.ones_like(sigma)
    # Sam's Meyer model
    elif site.site_model == 'sam' or site.site_model == 'Sam':
        rate = bm.promoter_curve_Meyer(site.k_min, sigma)
    # Max-min model according oparams measured with SIDD
    elif site.site_model == 'maxmin' or site.site_model == 'Maxmin':
        rate = bm.promoter_curve_opening_E_maxmin(site.k_min, site.k_max, sigma, *site.oparams)
    # Inverted max-min model (where it is positive supercoiling sensitive)
    elif site.site_model == 'maxmin_I' or site.site_model == 'Maxmin_I':
        rate = bm.promoter_curve_opening_E_maxmin_I(site.k_min, site.k_max, sigma, *site.oparams)
    # Similar to max-min but with the effective energy
    elif site.site_model == 'effE' or site.site_model == 'EffE':
        rate = bm.promoter_curve_opening_E(site.k_min, sigma, sigma0=0, *site.oparams)
    elif site.site_model == 'none' or site.site_model == 'None' or site.site_model is None:
        print('This site ', site.name, ' does not have a binding model')
        rate = site.k_min * np.ones_like(sigma)
    else:  # If there's no model, there's no binding
        print('This site model ', site.site_model, ' cannot be recognized')
        rate = site.k_min * np.ones_like(sigma)
    return rate, sigma


# It computes the activity curves of topoisomerase activity with a continuum model.
# Input: 1.- topoisomerase. Optionals= sigma_min, sigma_max, delta_sigma, dt
# Output: 1.- Supercoiling density removed by time-step dt (topo_curve)
def topoisomerase_activity_curves_continuum(topo, sigma_min=-.2, sigma_max=.1, delta_sigma=.001, dt=1):
    sigma = np.arange(sigma_min, sigma_max, delta_sigma)
    if topo.name == 'topoI':
        topo_curve = em.topo1_continuum(sigma, topo.concentration, topo.k_cat, dt)
    elif topo.name == 'gyrase':
        topo_curve = em.gyrase_continuum(sigma, topo.concentration, topo.k_cat, dt)
    else:
        print('Could not recognize name of topoisomerase')
        topo_curve = np.zeros_like(sigma)
    return topo_curve, sigma


# TODO: Write what this function does
def topoisomerase_activity_curves_stochastic(topo, sigma_min=-.2, sigma_max=.1, delta_sigma=.001, dt=1):
    sigma = np.arange(sigma_min, sigma_max, delta_sigma)
    if topo.name == 'topoI':
        topo_curve = bm.topoI_binding(topo.k_on, topo.concentration, sigma)
    elif topo.name == 'gyrase':
        topo_curve = bm.gyrase_binding(topo.k_on, topo.concentration, sigma)
    else:
        print('Could not recognize name of topoisomerase')
        topo_curve = np.zeros_like(sigma)
    return topo_curve, sigma


# This function returns a very simple rate, which is calculated as the sum of events (average) in a time interval.
# Inputs: Initiation_signal = binding signal (1 when a initiation event took place, 0 when initiation didn't happen), &
# time, which is an array composed of the times (in general we only care about the last and initial times.
# Outputs: rate = which is simple the sum of initiation events divided by the time interval.
def simple_initiation_rate(initiation_signal, t0, tf):
    rate = np.sum(initiation_signal) / (tf - t0)
    return rate


# Function that calculates rates from the steady state.
# The steady state is assumed to be covered by the intervals time[ta,tb].
# This steady state is calculated with the curve = log( sum of initiation events / time ).
# When this curve reaches a plateau, it is assumed that the steady state was reached.
# Inputs: sites_df (the binding signal is taken), the time array, and the ranges [ta, tb] which indicate the interval
# of the steady state.
# Outputs: 1.- curves, which is the curve is the logarithmic curve we use to determine the steady state.
# 2.- rates, which are calculated from the plateau of the curve. 3.- the corresponding site names...
def initiation_rates(sites_df, time, ta=0, tb=-1):
    sites_names = sites_df.drop_duplicates(subset='name')['name']
    curves = []
    rates = []
    names = []
    for name in sites_names:
        mask = sites_df['name'] == name
        signal = sites_df[mask]['binding'].to_numpy()
        curve, rate = calculate_steady_state_initiation_curve(signal, time, ta, tb)
        curves.append(curve)
        rates.append(rate)
        names.append(name)
    return curves, rates, names


# Similar than the previous function, except it calculates the initiation rates by site type, e.g. genes.
# TODO: plotea el sumini y logsumini para ver la diferencia, anotalo a tu cuaderno/latex, y luego comentalo
def initiation_rates_by_type(sites_df, my_type, time, ta=0, tb=-1):
    mask = sites_df['type'] == my_type
    my_df = sites_df[mask]
    sites_names = my_df.drop_duplicates(subset='name')['name']
    curves = []
    rates = []
    names = []
    for name in sites_names:
        mask = my_df['name'] == name
        signal = my_df[mask]['binding'].to_numpy()
        curve, rate = calculate_steady_state_initiation_curve(signal, time, ta, tb)
        curves.append(curve)
        rates.append(rate)
        names.append(name)
    return curves, rates, names


# This function is the one that actually calculates the steady state initiation curve.
# Inputs: 1.- y=binding signal. 2.- time array. 3.- intervals [ta,tb] in which the plateau of the steady state is
# calculated
# Outputs: 1.- curve = steady state curve. 2.- The corresponding rate calculated from the curve
def calculate_steady_state_initiation_curve(y, time, ta, tb):
    frames = len(y)
    sum_ini = np.zeros(frames)  # this is the sum of the initiation events
    log_sum_ini = np.zeros(frames)  # And the logarithm of these sum of events divided by the time.
    for k in range(1, frames):
        sum_ini[k] = np.sum(y[0:k])

    if np.sum(sum_ini) <= 0.0:
        rate = 0.0000000001
    else:
        log_sum_ini = np.log(sum_ini / time)
        plateau = np.mean(log_sum_ini[ta:tb])
        rate = np.exp(plateau)
    curve = log_sum_ini
    return curve, rate

# TODO: Relative initiation rate is next, what is that?

# TODO: Make the following functions:
#  In the plotting script, maybe I should load the circuit and the DFs - not this script
#  1.- Binding signals:
#   1.1.- Build signal - Give it a site #enzyme - return signal with 0s/1s (frames)
#   1.2.- Build signals by type - Give sites_df - Return array of list with signals and list with names
#   1.3.- Build signals of all sites_df - Give sites_df - Returns array of (n_sites, frames)
#  2.- Compute cross-correlation hyper matrix: give signals and time t0
#      return matrix(n_sites,n_sites,frames) - diagonal is auto-correlations and off-diagonal cross-correlations
#  3.- Sites binding curves:
#     3.1.- Binding curve: Give model and calculate supercoiling sensitive curve.
#  4.- Supercoiling at promoter? Nah, this might be already there, but in my plotting script I should do it.
#  5.- Topoisomerase activity curves (for continuum)?
#  6.- Rates:
#  6.1.- Initiation rate in interval [t1,t2]
#  6.2.- Relative initiation rate
#  6.3.- Elongation rate/unbinding rate
#  6.4.- mRNA copy number
#  6.5.- mRNA synthesis rate
#  6.6.- mRNA relative synthesis rate
#  7.- Superhelical distribution at site - give sigma at site, return distribution.
