import numpy as np


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
