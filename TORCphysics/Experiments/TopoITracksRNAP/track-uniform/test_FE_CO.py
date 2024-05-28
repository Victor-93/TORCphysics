import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from TORCphysics import Circuit
from TORCphysics import topo_calibration_tools as tct


# Description
# ---------------------------------------------------------
# I will try to compute the average Fold-Enrichment (FE) and the correlation between the RNAP and topo I
# Additionally, I will plot the Fold-Enrichment curve and kde for the RNAP.
# The intention of this script is to visually test the procedures performed in the calibration process.

# Inputs
# ---------------------------------------------------------
reference_topoI_file = '../noRNAP/reference_topoI.txt'  # This one has the kde
target_gene_name = 'reporter'  # This for calculating correlation
reference_path = '../noRNAP/'

# Simulation conditions - Even though we won't run the simulation again
# --------------------------------------------------------------
dt = 0.25
initial_time = 0
final_time = 1000
time = np.arange(initial_time, final_time + dt, dt)

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = '../sites.csv'
enzymes_filename = None
environment_filename = 'environment_calibration.csv'
output_prefix = 'noRNAP'
frames = len(time)
series = True
continuation = False

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'topoI': 'red', 'gyrase': 'cyan'}
colors_dict = {'RNAP': 'black', 'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}
ylabels = ['Density', 'Fold-Enrichment', 'Fold-Enrichment']

# Let's load the circuit so we have some info
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Get target site
target_gene = [site for site in my_circuit.site_list if site.name == target_gene_name][0]

# Define x-axes
x_system = tct.get_interpolated_x(1, my_circuit.size)
x_gene = tct.get_interpolated_x(target_gene.start, target_gene.end)

# We need to define arrays for comparison
FE_dict = {'topoI':None, 'gyrase':None}  # This ones are fold-enrichment curves
kde_gene = {'RNAP':None, 'topoI':None, 'gyrase':None}  # This ones the kdes extrapolated to the gene regions

# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True)
for p, name in enumerate(names):

    ax = axs[p]

    nbins = tct.calculate_number_nbins(my_circuit, name)

    # Load data
    x = np.loadtxt('position_'+name+'.txt')
    x = x[~np.isnan(x)]  # Just in case to remove nans

    # Calculate kde
    kde_x, kde_y = tct.calculate_KDE(x, nbins, scaled=True)

    # Different process for RNAP or topos
    if name == 'RNAP':
        y_interpolated = tct.get_interpolated_kde(kde_x, kde_y, x_gene)
        kde_gene[name] = y_interpolated

        # Plot histogram
        hist = sns.histplot(x, bins=nbins, ax=ax, color=colors_dict[name], label=name)

        # For plotting, just plot the interpolated x,y
        x = x_gene
        y = y_interpolated
    else:
        y_interpolated = tct.get_interpolated_kde(kde_x, kde_y, x_system)

        # We also need to interpolate the reference
        kde_ref = np.loadtxt(reference_path+'reference_'+name+'.txt')
        y_ref = tct.get_interpolated_kde(kde_ref[:,0] , kde_ref[:,1], x_system)

        # Compute FE
        FE_curve = y_interpolated / y_ref
        FE_dict[name] = FE_curve  # Store it in dict

        # And let's extrapolate to the gene region
        y_gene = tct.get_interpolated_kde(kde_x, kde_y, x_gene)
        kde_gene[name] = y_gene  # And save to the dict

        # For plotting, let's plot the Fold enrichment curves
        x = x_system
        y = FE_curve

    ax.plot(x, y, color=colors_dict[name], lw=lw)

    # Labels
    # ------------------------
    ax.set_ylabel(ylabels[p], fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(name, fontsize=title_size)
    #ax.set_ylim(0,np.max(kde_y))
    ax.set_xlim(0, my_circuit.size)
    ax.grid(True)

# Let's compute average fold-enrichment
avg_FE = FE_dict['topoI'].mean()

# And the correlation matrix
correlation_matrix = np.corrcoef(kde_gene['topoI'], kde_gene['RNAP'])
correlation = correlation_matrix[0, 1]

print('mean FE', avg_FE)
print('correlation', correlation)

plt.show()
#plt.savefig(out+'.png')
#plt.savefig(out+'.pdf')


