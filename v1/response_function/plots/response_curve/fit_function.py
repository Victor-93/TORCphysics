
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
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
    return a+b/(1 + np.exp( -(x - sigma_t)/epsilon ) )
#def opening_energy(x, a, b, c):
#    return  c*((x-a)/b)/(1 + np.abs( (x-a)/b ) )
    

#Initial conditions
#--------------------------------------------------------------------------------------------
f0 = [-0.004, 0.005, 5, 4]
path="../../filtering/"

width = 6
height = 3

bp0 = 99
bp1 = 113

coutput="fitting_Gop"
#Load and prepare sequence
#--------------------------------------------------------------------------------------------
#Info file to read the data
sequence_f = path+"sequence.txt"

#Let's read the file
f = open( sequence_f, "r")
sequence = [] #With characters
for x in f:
    sequence.append(x[0])

#seq_bp = sequence[bp]

#Load and prepare the heatmaps (P and G)
#--------------------------------------------------------------------------------------------
#Info file to read the data
info = np.loadtxt(path+"info.txt")

#I need a test file to allocate the data to the heatmap.
#I assume that there will always be a G_0.txt
test = np.loadtxt( path+"G_0.txt" )

nbp = len(test[:,0])      #number of bp
n_sigma = len(info[:,0])  #number of supercoiling densities

#Initialize arrays
G_map = np.zeros( ( n_sigma, nbp ) )
P_map = np.zeros( ( n_sigma, nbp ) )
sigma = np.zeros( n_sigma )

#And read my data
for i in  range(n_sigma):

    sigma[i] = info[i,1]   #supercoiling density

    #For  G
    #----------------------
    cdat = path + "G_" + str( int( info[i,0] ) ) + ".txt"

    test = np.loadtxt( cdat )

    G_map[i,:] = test[:,1]  #load energies at supercoiling density sigma[i]

    #For  P
    #----------------------
    cdat = path + "P_" + str( int( info[i,0] ) ) + ".txt"

    test = np.loadtxt( cdat )

    P_map[i,:] = test[:,1]  #load probabilities at supercoiling density sigma[i]


#Plot data
#--------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width,height), tight_layout=True)
#fig, axs = plt.subplots(3, figsize=(width,height*3), sharex=True, sharey=False, constrained_layout=True)

#Plot G
#------------------------------------------------------------------
ax = axs
x = sigma

#Calculate the average and standard deviation fo the region I'm interested
aG = np.zeros(n_sigma)
sG = np.zeros(n_sigma)

for i in range(n_sigma):
    aG[i] = np.mean(G_map[i,bp0:bp1] )
    sG[i] = np.std(G_map[i,bp0:bp1] )

popt, pcov = curve_fit(opening_energy, x, aG, p0=f0)
print(popt)
#print(pcov)
np.savetxt("opening_energy_parameters.txt", popt)

ax.plot(x,aG, color="black", lw=2, label="SIDD")
ax.fill_between( x, aG+sG, aG-sG, color="black", alpha=0.5)

x = np.arange(-0.3,0.3, 0.001)
my_fit = opening_energy( x, *popt )
ax.plot(x, my_fit, color="blue", label="fit")


ax.set_ylabel(r'Opening energy (kcal/mol)')
ax.set_xlabel(r'$\sigma$')
ax.grid(True)
ax.legend(loc="best")

plt.savefig(coutput+".pdf")
plt.savefig(coutput+".png")


