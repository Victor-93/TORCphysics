import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import pandas as pd
# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# We want to visualize the outputs of the genearch experiments

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
input_data = 'genearch_experiment.pkl'
file_out = 'genearch_experiment'
promoter_case = 'medium'

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
fig, axs = plt.subplots(2,2, figsize=(2*width, 2*height), tight_layout=True)#, sharex=True)

my_results_list = []
# Let's plot some basic information. Let's plot transcription rate as a function of system size
for n, results in enumerate(data): # Go through simulation results

    # Unpack data
    info = results['system_info']
    rate_left = results['rate_left']
    rate_right = results['rate_right']

    # Calculate size of the system
    system_size = (info['dist1'] + info['gene_length'] + info['intergene_length'] +
                   info['gene_length'] + info['dist2'])

    # We need to sort the data in a clever way.
    # We need a dataframe that contains system_size = Small, Large, Intermediate 1 and Intermediate 2
    # Gene orientation: Divergent, Tandem, Convergent
    # Rate: value
    # Gene: Left, Right

    for i, rates in enumerate([rate_left, rate_right]):
        for j, rate in enumerate(rates):
            my_results_dict = {'Rate': float(rate)}

            if i == 0:
                my_results_dict['Gene'] = 'Left'
            else:
                my_results_dict['Gene'] = 'Right'

            # Sort system size
            if info['dist1'] == 200 and info['dist2'] == 200:
                my_results_dict['system_size'] = 'Small'
            if info['dist1'] == 1000 and info['dist2'] == 1000:
                my_results_dict['system_size'] = 'Large'
            if info['dist1'] == 1000 and info['dist2'] == 200:
                my_results_dict['system_size'] = 'Intermediate 1'
            if info['dist1'] == 200 and info['dist2'] == 1000:
                my_results_dict['system_size'] = 'Intermediate 2'

            # Sort system configuration
            if info['g1_orientation'] == -1 and info['g2_orientation'] == 1:
                my_results_dict['Configuration'] = 'Divergent'
            if info['g1_orientation'] == info['g2_orientation']:
                my_results_dict['Configuration'] = 'Tandem'
            if info['g1_orientation'] == 1 and info['g2_orientation'] == -1:
                my_results_dict['Configuration'] = 'Convergent'

            my_results_list.append(my_results_dict)

results_df = pd.DataFrame(my_results_list)
    # Plot the data
    #axs.errorbar(system_size, np.mean(rate_left), yerr=np.std(rate_left), color='red')
    #axs.errorbar(system_size, np.mean(rate_right), yerr=np.std(rate_right), color='blue')

# TODO: Change the conditions here for selecting the systems of interest
axs[0,0].set_title('Large domain', fontsize=title_size)
sns.barplot(results_df.loc[results_df['system_size'] == "Large"], x="Configuration", y="Rate", hue="Gene", ax=axs[0,0])
#sns.barplot(results_df, x="Configuration", y="Rate", hue="Gene", ax=axs[0,0])


# Sort labels
for ax in axs.flatten():
    ax.grid(True)
    ax.set_xlabel('Gene configuration', fontsize=xlabel_size)
    ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)

# TODO: You can save in png format by uncommenting the png line
#plt.savefig(file_out+'.png')
#plt.savefig(file_out+'.pdf')

plt.show()