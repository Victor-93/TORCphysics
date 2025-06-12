from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import re
import sys


#Description
#--------------------------------------------------------------------------------------------
#I am going to calculate the average G of the openning energy by averaging
#the -10 region, the discriminator and the +3 bp... Note that the bubble takes a little bit
#of both the -10 (the machine learning paper says that half of the -10 region), but here I'll
#consider all of it as I don't think it'll change much.

#Functions
#--------------------------------------------------------------------------------------------

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


#def opening_energy(x, a, b, c):
#    return  c*((x-a)/b)/(1 + np.abs( (x-a)/b ) )


#Initial conditions
#--------------------------------------------------------------------------------------------
case = 'medium'

f0 = [0.09, 2.00, -.1, .007]
path = "../../filtering/"

width = 6
height = 3

sigma0 = 10#150  # should be -15

title = "gene architecture - medium promoter"

coutput = "fitting_Gop"

alpha = 0.15
lw = 1

x0 = np.arange(-0.2, 0., 0.001)

#Gene colours
weak_c = 'green'
medium_c = 'blue'
strong_c = 'red'
#Sequences
#--------------------------------------------------------------------------------------------
# Promoter sequences
weak_seq ='AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGATAGCGG'
medium_seq = 'TTGACATCAGGAAAATTTTTCTGCATAATTATTTCATATCAC'
strong_seq = 'TTGACATCGCATCTTTTTGTACCTATAATGTGTGGATAGAGT'

#Load and prepare sequence
#--------------------------------------------------------------------------------------------
#Info file to read the data
sequence_f = path + "sequence.txt"

#Let's read the file
f = open(sequence_f, "r")
sequence = []  #With characters
for x in f:
    sequence.append(x[0])

#Locate regions of interest
#--------------------------------------------------------------------------------------------
plasmid_sequence = ''.join(sequence)

#Locate regions of interest
#--------------------------------------------------------------------------------------------
plasmid_sequence = ''.join(sequence)

x_Pweak =  plasmid_sequence.find(weak_seq)
n_Pweak = len(weak_seq)

x_Pmedium = plasmid_sequence.find(medium_seq)
n_Pmedium = len(medium_seq)

x_Pstrong = plasmid_sequence.find(strong_seq)
n_Pstrong = len(strong_seq)

#Load and prepare the heatmaps (P and G)
#--------------------------------------------------------------------------------------------
#Info file to read the data
info = np.loadtxt(path + "P_M.txt")  #Let's treat this as info

#I need a test file to allocate the data to the heatmap.
#I assume that there will always be a G_0.txt
test = np.loadtxt(path + "G_M_0.txt")

nbp = len(test[:, 0])  #number of bp
n_sigma = len(info[:, 0])  #number of supercoiling densities

#Initialize arrays
G_map = np.zeros((n_sigma, nbp))
P_map = np.zeros((n_sigma, nbp))
sigma = np.zeros(n_sigma)

#And read my data
for i in range(n_sigma):
    sigma[i] = info[i, 1]  #supercoiling density

    #For  G
    #----------------------
    #cdat = path + "G_M_" + str(int(info[i, 0])) + ".txt"
    cdat = path + "G_M_" + str(i) + ".txt"

    test = np.loadtxt(cdat)

    G_map[i, :] = test[:, 1]  #load energies at supercoiling density sigma[i]

G_map[G_map < 0] = 0  # Filter values lower than 0. This is a mistake! Why did it happen?

sigma = info[:, 0]
#Plot data
#--------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
#fig, axs = plt.subplots(3, figsize=(width,height*3), sharex=True, sharey=False, constrained_layout=True)

fig.suptitle(title)

#Plot G
#------------------------------------------------------------------
ax = axs
x = sigma

if case == 'weak':
    bp0 = x_Pweak
    bp1 = x_Pweak + n_Pweak
    mcolor = weak_c
    cparams = "opening_energy_parameters_" + case + '.txt'
    label = case

if case == 'medium':
    bp0 = x_Pmedium
    bp1 = x_Pmedium + n_Pmedium
    mcolor = medium_c
    cparams = "opening_energy_parameters_" + case + '.txt'
    label = case

if case == 'strong':
    bp0 = x_Pstrong
    bp1 = x_Pstrong + n_Pstrong
    mcolor = strong_c
    cparams = "opening_energy_parameters_" + case + '.txt'
    label = case

#Calculate the average and standard deviation fo the region I'm interested
aG = np.zeros(n_sigma)
sG = np.zeros(n_sigma)

for i in range(n_sigma):
    aG[i] = np.mean(G_map[i, bp0:bp1])
    sG[i] = np.std(G_map[i, bp0:bp1])

print(aG[:])
print(x[:])
popt, pcov = curve_fit(opening_energy, x[sigma0:], aG[sigma0:], p0=f0)
np.savetxt(cparams, popt)

ax.plot(x, aG, color=mcolor, lw=lw, label=label)
ax.fill_between(x, aG + sG, aG - sG, color=mcolor, alpha=alpha)

my_fit = opening_energy( x0, *popt )
ax.plot(x0, my_fit, "--", color=mcolor, lw=lw, label="fit")

# Ticks stuff
#-----------------------------------------------------------------------------
ax.set_ylabel(r'Opening energy (kcal/mol)')
ax.set_xlabel(r'$\sigma$')
ax.grid(True)
ax.legend(loc="best")
ax.set_ylim(-0.5, 13)

plt.savefig(coutput + ".pdf")
plt.savefig(coutput + ".png")
