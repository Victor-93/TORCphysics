from TORCphysics import Circuit
from TORCphysics import visualization as vs
import matplotlib.pyplot as plt
import numpy as np

# Simulation conditions
circuit_filename = 'circuit.csv' # Linear version
#circuit_filename = 'circuit_circular.csv' # Circular version
sites_filename = 'single_gene_site.csv'
enzymes_filename = None  # This input is optional, no enzymes_file means that nothing is bound
#environment_filename = 'environment_torcphys.csv'  # In this one, topoisomerases can bind (more realistic but slower).
                                                   # It similar than environment_3.csv but uses the topoI and gyrase params csv ile
environment_filename = 'environment_4.csv'  # Topoisomerases act continuously in this one (they dont bind so its faster)
output_prefix = 'single_gene_continuum-topos'
frames = 5000
#frames = 200 # Just a few for the video - specially if topoisomerase can bind as we have more objects binding to the DNA
series = True
continuation = False
dt = 1.0
do_animation = False
animation_format = '.gif' # or '.mp4'
#animation_format = '.mp4' # or '.mp4'

colors_dict = {'single': 'purple', 'left': 'blue', 'right': 'red'}

# Initialize circuit with the initial conditions
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
# Print general information
#my_circuit.print_general_information()

# And run simulation
my_circuit.run()

# Do the animation
if do_animation:
    # Produce animation either as a gif or a mp4 - If producing a movie, maybe reduce the number of frames since it could take a while.
    vs.create_animation_linear_artist(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, my_circuit.environmental_df,
                               my_circuit.name, out_format=animation_format,
                               site_type='gene', site_colours=colors_dict, draw_containers=False,fps=30)

# *********************************************
# If using pycharm, with the debug option you can check the elements on the my_circuit object and dataframes


# Let's do some plots
fig, axs = plt.subplots(2, figsize=(7, 3*2), tight_layout=True)

fig.suptitle('Single short gene')

# Plot the signal profiles, just in case
vs.plot_signal_profiles(my_circuit=my_circuit, sites_df=my_circuit.sites_df, site_type='gene', colors=colors_dict, axs=axs[0])

# Let's plot the transcription rate. We can obtain it by dividing the number of mRNA produced by the time
gene_name = 'single'
gene_mask = my_circuit.environmental_df['name'] == gene_name # This gives the positions for which the condition is met
mRNA = my_circuit.environmental_df.loc[gene_mask] # This is the filtered dataframe

axs[1].plot(mRNA['time'], mRNA['concentration']/mRNA['time']) # Plot the number of accumulated mRNAs/time (mRNA should degrade in real life, but we don't include it in the model)

average_rate = np.mean(mRNA['concentration']/mRNA['time']) # Very vaguely, we can calculate the transcription rate as an average. Ideally we should avoid the first spikes....
average_rate_s = "{:.2e}".format(average_rate)  # Let's save it as a string

# We can print the average rate to the screen or the plot
print('Average transcription rate:', average_rate_s)
axs[1].text(0.2,0.8, r'$k=$'+average_rate_s, ha='center', transform=axs[1].transAxes, fontsize=12) # This is to print it in the plot


# Make the graph nice
axs[1].set_xlabel('Time')
axs[1].set_ylabel('mRNA/time')
axs[1].grid(True)

#plt.savefig('single_gene.png') # We can save the figure as well
plt.show()
