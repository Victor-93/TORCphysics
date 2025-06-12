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
promoters = ['Weak', 'Medium', 'Strong']
paths = ['weak/filtering/', 'medium/filtering/', 'strong/filtering/']
promoter_colors = ['green', 'blue', 'red']
title = "gene architecture - strong promoter"

circle = False #if circle
w=16 #window size for calculating gc content

#Plotting variables
coutput = "energy_profiles"
#P_colour = 'Blues_r'
my_colour = "blue"
GC_colour = 'black'
aspectr  = 1.5#.5 #40
width = 5
height = 2.5#3.75
lw = 4
font_size = 12
xlabel_size = 14
title_size = 16
ms=6

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

y_genes = [5, 5]  # For heatmaps
#y_genes = [195, 195]  # For heatmaps
y_genes2 = [99, 99] #The one used for the GC contenty
g_lw = 6 #Genes line width

vars = ['G_M']#, 'Prob_M']#, 'G_Z', 'G_C']
colour_maps = ['magma_r']#, 'inferno']#, 'cool', 'viridis']
labels = ['Melting Energy (kcal/mol)']#, 'Melting probability']#, 'Cruciform probability', 'Z-DNA probability' ]

#Sequences
#--------------------------------------------------------------------------------------------
# Promoter sequences
weak_seq ='AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGATAGCGG'
medium_seq = 'TTGACATCAGGAAAATTTTTCTGCATAATTATTTCATATCAC'
strong_seq = 'TTGACATCGCATCTTTTTGTACCTATAATGTGTGGATAGAGT'

promoter_sequences = [weak_seq, medium_seq, strong_seq]
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

def plot_genes(ax, x_P,n_P, mycolor, ygene):
    ax.plot( [x_P, x_P+n_P], ygene, lw=g_lw, color=mycolor)

#PROCESS - Let's plot as we process the data
#--------------------------------------------------------------------------------------------
fig, axs = plt.subplots(len(promoters), figsize=(width,height*len(promoters)), sharex=True, sharey=False, tight_layout=True)

for p, pcase in enumerate(promoters):
    path = paths[p]
    promoter_sequence = promoter_sequences[p]
    ax = axs[p]
    ax.set_title(promoters[p] + ' Promoter', color=promoter_colors[p], fontsize=title_size)

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
    genome_sequence = ''.join(sequence)

    x_P = genome_sequence.find(promoter_sequence)
    n_P = len(promoter_sequence)

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

    #Plot vars
    #------------------------------------------------------------------
    for i, myvar in enumerate(vars):
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

        plot_genes(ax, x_P, n_P, promoter_colors[p], y_genes)

        ax.set_yticks( y_pos )
        ax.set_yticklabels( y_sigma)
        ax.set_ylabel(r'Superhelical Density', fontsize=xlabel_size)
        ax.grid(True)
axs[2].set_xlabel('Position (bp)', fontsize=xlabel_size)
plt.subplots_adjust(hspace=.1)
#plt.savefig(coutput+".pdf")
plt.savefig(coutput+".png")
plt.show()



