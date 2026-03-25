import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import pickle
from TORCphysics import visualization as vs
from TORCphysics import analysis as an
# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# We want to visualize the outputs of the genearch experiments

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
input_data = 'genearch_experiment.pkl'
file_out = 'genearch_experiment'
promoter_case = 'weak'

# -----------------------------------
# FIGURE Params
# -----------------------------------
width = 7
height = 3.5
lw = 3
font_size = 12
xlabel_size = 16
title_size = 18
nbins=30

colors = ['black', 'red', 'blue', 'green', 'yellow', 'purple', 'gray']

# ----------------------------------------------------------------------------------------------------------------------
# Load simulation outputs
# ----------------------------------------------------------------------------------------------------------------------
with open(input_data, 'rb') as file:
    # Load the pickled data
    data = pickle.load(file)  # This data contains information regarding the experiments
                              # It is a list of dictionaries. Each list entry represents the results from the condiitons

# USE THE DEBUGGING TOOL HERE TO EXPLORE THE CONTENTS

# ----------------------------------------------------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)#, sharex=True)

# Let's plot some basic information. Let's plot transcription rate as a function of system size
for n, results in enumerate(data): # Go through simulation results

    # Unpack data
    info = results['system_info']
    rate_left = results['rate_left']
    rate_right = results['rate_right']

    # Calculate size of the system
    system_size = (info['dist1'] + info['gene_length'] + info['intergene_length'] +
                   info['gene_length'] + info['dist2'])

    # Plot the data
    axs.errorbar(system_size, np.mean(rate_left), yerr=np.std(rate_left), color='red')
    axs.errorbar(system_size, np.mean(rate_right), yerr=np.std(rate_right), color='blue')

# Sort labels
axs.set_ylabel('Transcription rate (transcript/s)')
axs.grid(True)
axs.set_xlabel('System size (bp)', fontsize=xlabel_size)

#plt.savefig(file_out+'.png')
#plt.savefig(file_out+'.pdf')

plt.show()