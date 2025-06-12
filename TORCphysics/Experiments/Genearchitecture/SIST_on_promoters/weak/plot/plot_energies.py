import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from csv import reader
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import re
import sys

#Description
#--------------------------------------------------------------------------------------------
#I need to plot the sequence and heatmaps of P and G

#Initial conditions
#--------------------------------------------------------------------------------------------
path="../filtering/"

title = "gene architecture - weak promoter"

circle = False #if circle
w=16 #window size for calculating gc content

#Plotting variables
coutput = "profiles"
#P_colour = 'Blues_r'
my_colour = "blue"
GC_colour = 'black'
aspectr  = 2 #40
width = 10#8
height = 4#3.5

#Gene colours
weak_c = 'green'
medium_c = 'blue'
strong_c = 'red'

#Drawing genome colours
ah = 115         #height
awidth=40.0      #width of arrow
gene_color='red' #gene color
#y_pos = [0, 50, 100, 150, 200]
#y_sigma = [-0.20, -0.15, -0.10, -0.15, 0.0]
y_pos = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
y_sigma = [0.0, -0.02, -0.04, -0.06, -0.08, -0.1, -0.12, -0.14, -0.16, -0.18]
#y_sigma = [0.0, -0.02, -0.04, -0.06, -0.08]
#y_sigma = [-0.08, -0.06, -0.04, -0.02, 0.0]

y_genes = [195, 195]  # For heatmaps
y_genes2 = [99, 99] #The one used for the GC contenty
g_lw = 6 #Genes line width

vars = ['G_M', 'Prob_M']#, 'G_Z', 'G_C']
colour_maps = ['magma', 'inferno']#, 'cool', 'viridis']
labels = ['Melting energy', 'Melting probability']#, 'Cruciform probability', 'Z-DNA probability' ]

#Sequences
#--------------------------------------------------------------------------------------------
# Promoter sequences
weak_seq ='AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGATAGCGG'
medium_seq = 'TTGACATCAGGAAAATTTTTCTGCATAATTATTTCATATCAC'
strong_seq = 'TTGACATCGCATCTTTTTGTACCTATAATGTGTGGATAGAGT'

#Functions
#--------------------------------------------------------------------------------------------
def gc_content(seq):
    #How much GC a sequence has.
    return float( seq.count('C') + seq.count('G') )/len(seq) * 100

def gc_content_window( seq, window=50 ):
    #How much content in a region
    #window = 30 by default because I assume that's the size of the RNAPolymerase
    content = []
    for i in range(0, len(seq) - window +1, 1): #This will scan the sequence
        subseq = seq[i:i + window]
        content.append( gc_content(subseq) )
    return content

def plot_genes(ax, ygene):
    ax.plot( [x_Pweak, x_Pweak+n_Pweak], ygene, lw=g_lw, color=weak_c)
    ax.plot( [x_Pmedium, x_Pmedium+n_Pmedium], ygene, lw=g_lw, color=medium_c)
    ax.plot( [x_Pstrong, x_Pstrong+n_Pstrong], ygene, lw=g_lw, color=strong_c)

#Load and prepare sequence
#--------------------------------------------------------------------------------------------
#Info file to read the data
sequence_f = path+"sequence.txt"

#Let's read the file
f = open( sequence_f, "r")
sequence = [] #With characters
for x in f:
    sequence.append(x[0])

nbp = len(sequence) #This is the real sequence, but for calculating the whole
                    #GC content we need to add elements so it can complete the circle

sequence_c = sequence

if circle:
    a=0
    for i in range(w-1):   #100 is the size of the window
        sequence_c.append(sequence[a])
        a=a +1

#Let's calculate the content
content = gc_content_window(sequence_c)
content = np.roll(content,int(w/2) ) #It needs to move half of the window so it representes the sourroundins.

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
info = np.loadtxt(path+"P_M.txt")  #Let's treat this as info

#I need a test file to allocate the data to the heatmap.
#I assume that there will always be a G_0.txt
test = np.loadtxt( path+vars[0]+"_0.txt" )

nbp = len(test[:,0])      #number of bp
n_sigma = len(info[:,0])  #number of supercoiling densities
print(nbp, n_sigma)

#Initialize arrays
bp_ax = test[:,0]
sigma = np.zeros( n_sigma )
my_map = np.zeros( ( n_sigma, nbp ) )

my_maps = []
my_maps = np.zeros( ( len(vars), n_sigma, nbp ) )
# Lets load the data
s = -1
for myvar in vars:
    s=s+1

    #my_map[:,:] = 0.0

    #And read my data
    for i in  range(n_sigma):

        sigma[i] = info[i,1]   #supercoiling density

        #Load
        #----------------------
        #cdat = path + myvar + "_" + str( int( info[i,0] ) ) + ".txt"
        cdat = path + myvar + "_" + str( i ) + ".txt"

        #print(cdat)

        test = np.loadtxt( cdat )

        my_maps[s,i,:] = test[:,1]
        my_map[i,:] = test[:,1]  #load energies/probabilities at supercoiling density sigma[i]
        #if i == 14 and myvar == 'C':
        #    print(test[3114,1], test[3114,0])
        #    sys.exit()

    #my_maps.append( my_map )

ss = 6 # half the supercoiling

#Plot data
#--------------------------------------------------------------------------------------------
#fig, axs = plt.subplots(2, figsize=(width,height*2), sharex=True, sharey=False)#, tight_layout=True)
fig, axs = plt.subplots(1+len(vars), figsize=(width,height*(1+len(vars))), sharex=True, sharey=False, constrained_layout=True)

fig.suptitle(title)

#Plot GC content
#------------------------------------------------------------------
ax = axs[0]
x = np.arange(1,len(content)+1)
y = content
ax.fill_between(x, y,color=GC_colour, alpha=0.5)
ax.set_ylabel(r'GC content (%)')
ax.set_xlabel('bp position')

plot_genes(ax, y_genes2)

ax.set_ylim(0,100)
ax.grid(True)


#Plot vars
#------------------------------------------------------------------
for i, myvar in enumerate(vars):
    ax = axs[i+1]
    my_map = my_maps[i,:,:]

    if 'G_' in myvar:
        im  = ax.imshow( my_map, cmap= colour_maps[i], interpolation= 'nearest', aspect=aspectr, vmin=0, vmax=12 )
    else:
        im  = ax.imshow( my_map, cmap= colour_maps[i], interpolation= 'nearest', aspect=aspectr)

    axins = inset_axes(ax,
                   width="2.5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
    fig.colorbar(im, cax=axins, label=labels[i])# ticks=[1, 2, 3])

    plot_genes(ax, y_genes)

    ax.set_yticks( y_pos )
    ax.set_yticklabels( y_sigma)
    ax.set_ylabel(r'$\sigma$')
    ax.set_xlabel('bp position')
    ax.grid(True)


plt.subplots_adjust(hspace=.1)

plt.savefig(coutput+".pdf")
plt.savefig(coutput+".png")



