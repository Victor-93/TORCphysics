import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#Description
#-----------------------------------------------------------------------------------------------------------------------
# I want to create a supplementary figure with the three promoter responses and their

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
path = '/home/victor/Documents/Research/TORCPhyspaper/gene_architecture/SIST_on_promoters/'
promoter_cases = ['weak', 'medium', 'strong']

SIST_paths = []
for pcase in promoter_cases:
    SIST_paths.append(path+pcase+'/filtering/')

promoter_sequences = [
    'AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGATAGCGG',
    'TTGACATCAGGAAAATTTTTCTGCATAATTATTTCATATCAC',
    'TTGACATCGCATCTTTTTGTACCTATAATGTGTGGATAGAGT'
]

sigma0 = 10#150  # should be -15
f0 = [0.09, 2.00, -.1, .007]

#Figure params
#-----------------------------------------------------------------------------------------------------------------------
width = 6.5
height = 3#3.75
lw = 2
font_size = 12
xlabel_size = 14
title_size = 16
alpha=0.25

# line styles
model_ls = '-o'
exp_ls = '--o'
titles = ['Weak Promoter', 'Medium Promoter', 'Strong Promoter']

colors = ['green', 'blue', 'red']

#Functions
#-----------------------------------------------------------------------------------------------------------------------

#The opening energy which follows a sigmoidal curve.
#x = supercoiling density
#a = thresshold
#b = width of the crossover
#c = thermal constant
#d = I'm just curious

#def opening_energy(x, a, b, c):
#    return c/(1 + np.exp( -(x - a)/b ) )

def opening_energy(x, a, b, sigma_t, epsilon):
    #    return a + b * (x-c)/(1 + np.abs(x-c) )
    return a + b / (1 + np.exp(-(x - sigma_t) / epsilon))

#Process
#-----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True, sharex=True)

x0 = np.arange(-0.2, 0., 0.001)
promoter_label = ['Weak Promoter', 'Medium Promoter', 'Strong Promoter']
outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

for k, pcase in enumerate(promoter_cases):

    ax = axs[k]
    path = SIST_paths[k]
    ax.set_title(promoter_label[k], fontsize=title_size, color=colors[k])

    # Load and prepare sequence
    # --------------------------------------------------------------------------------------------
    # Info file to read the data
    sequence_f = path + "sequence.txt"  # This is the total sequence of the system

    my_seq = promoter_sequences[k]  # This is the sequence of my promoter

    # Let's read the file
    f = open(sequence_f, "r")
    sequence = []  # With characters
    for x in f:
        sequence.append(x[0])

    # Locate regions of interest
    # --------------------------------------------------------------------------------------------
    plasmid_sequence = ''.join(sequence)

    x_P = plasmid_sequence.find(my_seq)
    n_P = len(my_seq)

    # Load and prepare the heatmaps (P and G)
    # --------------------------------------------------------------------------------------------
    # Info file to read the data
    info = np.loadtxt(path + "P_M.txt")  # Let's treat this as info

    # I need a test file to allocate the data to the heatmap.
    # I assume that there will always be a G_0.txt
    test = np.loadtxt(path + "G_M_0.txt")

    nbp = len(test[:, 0])  # number of bp
    n_sigma = len(info[:, 0])  # number of supercoiling densities

    # Initialize arrays
    G_map = np.zeros((n_sigma, nbp))
    P_map = np.zeros((n_sigma, nbp))
    sigma = np.zeros(n_sigma)

    # And read my data
    for i in range(n_sigma):
        sigma[i] = info[i, 1]  # supercoiling density

        # For  G
        # ----------------------
        # cdat = path + "G_M_" + str(int(info[i, 0])) + ".txt"
        cdat = path + "G_M_" + str(i) + ".txt"

        test = np.loadtxt(cdat)

        G_map[i, :] = test[:, 1]  # load energies at supercoiling density sigma[i]

    G_map[G_map < 0] = 0  # Filter values lower than 0. This is a mistake! Why did it happen?

    sigma = info[:, 0]

    # Plot G
    # ------------------------------------------------------------------
    bp0 = x_P
    bp1 = x_P + n_P
    mcolor = colors[k]
    x=sigma

    # Calculate the average and standard deviation fo the region I'm interested
    aG = np.zeros(n_sigma)
    sG = np.zeros(n_sigma)

    for i in range(n_sigma):
        aG[i] = np.mean(G_map[i, bp0:bp1])
        sG[i] = np.std(G_map[i, bp0:bp1])

    print(aG[:])
    print(x[:])
    popt, pcov = curve_fit(opening_energy, x[sigma0:], aG[sigma0:], p0=f0)

    ax.plot(x, aG, color=mcolor, lw=lw, label=r'$\bar{G}$')
    ax.fill_between(x, aG + sG, aG - sG, color=mcolor, alpha=alpha)

    my_fit = opening_energy(x0, *popt)
    ax.plot(x0, my_fit, "-", color='black', lw=lw, label=r"$U_{melt}$")

    # Add label outside the plot
    ax.text(-0.15, 1.15, outside_label[k], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

    # Ticks stuff
    # -----------------------------------------------------------------------------
    ax.set_ylabel(r'Melting Energy (kcal/mol)', fontsize=xlabel_size)
    ax.grid(True)
    ax.set_ylim(-0.5, 13)
axs[0].legend(loc="best", fontsize=font_size*1.2)
axs[2].set_xlabel(r'Superhelical Density', fontsize=xlabel_size)

plt.savefig('promoter_responses_fit.png')
#plt.savefig('promoter_responses_fit.pdf')
plt.show()
