import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import pandas as pd
import params
import statistical_model as sm
#from scipy import signal
import sys


#--------------------------------------------------------------------------
#DESCRIPTION
#--------------------------------------------------------------------------
#Ok, not it is time to plot averages

#params
#---------------------------------------------------------
width = 6
height = 2

n_sim = 10 #Number of simulations

cpath="results/" #Path of results

df = pd.read_csv( 'input.csv', sep='\t')

genome_df = df[ df[ 'type' ].isin(['genome']) ]          #dataframe of genome (to get size)
genome_df = genome_df.reset_index(drop=True)

genes_df  = df[ df[ 'type' ].isin(['gene']) ]            #dataframe of genes (to get orientations...)
genes_df = genes_df.reset_index(drop=True)

n_genes = len(genes_df['start']) #number of genes

o_df  = df[ df[ 'type' ].isin(['NAP']) ]   #dataframe of objects (NAPs )
n_NAPs = len(o_df['start'])           #number of objects bound to the DNA (NAPs)

#Let's load the output
output_array = np.load(cpath+'output_1.npy') #just one to get the number of frames
#output_array = np.load('output.npy')

frames = len(output_array[:,0,0] )

time = np.arange( 0, frames*params.dt, params.dt)

#Read global sigma
#globalsigma = np.loadtxt('results/global_sigma_5.txt')
#globalsigma = globalsigma[:,0]


colors = []
colors.append("green")
colors.append("brown")
colors.append("blue")
colors.append("magenta")
colors.append("yellow")
colors.append("purple")
colors.append("cyan")
colors.append("black")

#---------------------------------------------------------
#---------------------------------------------------------
#FIGURE
#---------------------------------------------------------
#---------------------------------------------------------

fig, axs = plt.subplots(2, figsize=(width*1.5, height*2.), tight_layout=True)#, sharex=True)

#We will plot 2 things:
#   1- system representation
#   2- average cross correlation

#---------------------------------------------------------
#Draw genes
#---------------------------------------------------------
ax = axs[0]
h=1.5
dh=0.5
gene_lw   = 5
object_text = 10 #NAPs, genes, etc...
DNA_lw    = 2
DNA_colour = 'black'


g0 = genome_df.iloc[0]['start']
g1 = genome_df.iloc[0]['end']
gx = np.array( [g0, g1] )
gy = np.array( [h, h] )
ax.set_xlim(g0-100,g1+100)
ax.set_ylim(0,2)
#ax.grid(axis='x', zorder=1)
ax.tick_params( labelleft=False, bottom=False,top=False)
ax.plot( gx,gy, lw=DNA_lw, color=DNA_colour, zorder=2)

for i in range(n_genes):
    dx = genes_df.iloc[i]['end'] - genes_df.iloc[i]['start']
    x0 = genes_df.iloc[i]['start']
    x1 = genes_df.iloc[i]['end']
    name = genes_df.iloc[i]['name']
    rate = genes_df.iloc[i]['rate']
    #arrow = mpatches.FancyArrowPatch( (x0,h-.01), (x1,h-.01), mutation_scale=25, color=gene_colour, zorder=3, lw=gene_lw)
    arrow = mpatches.FancyArrowPatch( (x0,h), (x1,h), mutation_scale=20, color=colors[i], zorder=3, lw=gene_lw)
    ax.add_patch(arrow)

    if x0 < x1:
        a = x0 + abs(dx/2)
    else:
        a = x1 + abs(dx/2)
    ax.text(a,h-dh, name, fontsize=object_text) #name
    ax.text(a,h-1.5*dh, r'$k_0=$'+str(rate), fontsize=object_text) #rate

#Draw naps
print(n_NAPs)
for i in range(n_NAPs):
    x0 = o_df.iloc[i]['start']
    print(x0)
    rect = mpatches.Rectangle( (x0,h-h/10), 21*2.5, h/5, color='black', lw=1) 
    #linewidth=1, edgecolor='r', facecolor='none')

    # Add the patch to the Axes
    ax.add_patch(rect)

ax.set_xlabel( 'bp position')

#---------------------------------------------------------
#PLOT CROSS CORRELATIONS 
#----------------------------------------------------------

#Let's First create our signals
#--------------------------
signal = np.zeros( ( frames, n_genes, n_sim ) )

for n in range(n_sim):
    output_array = np.load(cpath+'output_'+str(n+1)+'.npy') 
    for i in range(n_genes):
        for k in range(frames):
            if output_array[k,i,0] > 0: #If there is one RNAP transcribing
                signal[k,i,n] = 1

#This t0 is the time at which the cross-correlation takes place.
t0 = 1000 #600#10 mins of equilibration?
N = len(signal[t0:,0,0])

#This lag is the x axis of the cross correlation plot
lag = np.arange( -N*params.dt*.5, N*params.dt*.5, params.dt)

#Let's compute correlations
#--------------------------
ax = axs[1]
acorr = np.zeros( (N, n_genes, n_sim) )  #autocorrelations
crossc = np.zeros( (N, n_genes, n_sim) ) #cross-correlations with tetA
ncrossc = np.zeros( (N, n_genes, n_sim) ) #cross-correlations with tetA normalized

#Let's first compute auto correlations
for n in range(n_sim):
    for i in range(n_genes):
        acorr[:,i,n] = np.correlate( signal[t0:,i,0], signal[t0:,i,0],"same" ) #this is the autocorrelation

#Now the cross-correlations with tetA gene (in this case is index 0) 
for n in range(n_sim):
    for i in range(n_genes):
        crossc[:,i,n] = np.correlate( signal[t0:,0,n], signal[t0:,i,n],"same" ) #It is not normalized

        #And normalize
        if np.max(acorr[:,i,n]) > 0:
            ncrossc[:,i,n] = crossc[:,i,n]/np.sqrt( np.max(acorr[:,0,n])*np.max(acorr[:,i,n]) ) # This normalizes them
        else:
            ncrossc[:,i,n] = crossc[:,i,n]

#Now let's compute the averages for each timelag T
avg_ncrossc = np.zeros( (N,n_genes) )
std_ncrossc = np.zeros( (N,n_genes) )
for i in range(n_genes):
    for T in range(N):
        avg_ncrossc[T,i] = np.mean( ncrossc[T,i,:] )
        std_ncrossc[T,i] = np.std( ncrossc[T,i,:] )

#and PLOT THEM
for i in range(n_genes):
    if genes_df.iloc[i]['name'] == 'CDS':
        continue
    if genes_df.iloc[i]['name'] == 'tetA':
        continue
    a = avg_ncrossc[:,i]
    b = std_ncrossc[:,i]
    print( "max cross-correlation in gene", genes_df.iloc[i]['name'], " at (corr,timelag):", np.max(a), lag[np.argmax(a)] )
    ax.plot( lag, a, color=colors[i] )
    ax.fill_between( lag, a+b, a-b, color=colors[i], alpha=0.2   )
#ax.legend(loc='best')
ax.grid(True)
ax.set_xlim(-750,750)
ax.set_ylabel( 'Cross-corr w tetA')
ax.set_xlabel( 'time lag (seconds)')

plt.savefig('avg-cross-correlation.pdf')
plt.savefig('avg-cross-correlation.png')

