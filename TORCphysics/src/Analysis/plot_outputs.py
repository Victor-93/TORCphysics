import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import params

#from scipy import signal


#--------------------------------------------------------------------------
#DESCRIPTION
#--------------------------------------------------------------------------
#I'm going to plot the circle, and the outputs (transcripts).
# * Translation is assumed to start as soon as elongation starts and until it finishes.
# * It is also assumed that the proportion of fluorescent proteins created is proportional
#   to the number of transcripts created
#That's why we approximate the itensities (outputs) to the number of RNAPs transcribing the genes
#in each population

#params
#---------------------------------------------------------
width = 6
height = 3

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

time = np.arange(0, frames * params.dt, params.dt)

#Better use these ones for the colors of genes...
#Sort them according the input file...
colors = []
colors.append("yellow")
colors.append("black")
colors.append("blue")
colors.append("red")
colors.append("magenta")
colors.append("green")
colors.append("cyan")
colors.append("black")

NAP_ms = 8
ori_ms = 15
ori_s  = 'p'
NAP_s  = "s"
ori_c  = "yellow"   
NAP_c  = "gray"

#---------------------------------------------------------
#---------------------------------------------------------
#FIGURE
#---------------------------------------------------------
#---------------------------------------------------------

fig, ax = plt.subplots(1, figsize=(width, height), tight_layout=True)

#---------------------------------------------------------
#PLOT SUPERCOILING AND SIGNALS
#----------------------------------------------------------

#Let's First create our signals (which is not 
#--------------------------
signal = np.zeros( ( frames, n_genes ) )
sigma = np.zeros( ( frames, 2 ) ) #0 = avg, 1 = std - global sigma
s = np.zeros( ( frames, n_sim ) ) #this will help us calculate the average sigma

for n in range(n_sim):
    output_array = np.load(cpath+'output_'+str(n+1)+'.npy')
    globalsigma = np.loadtxt(cpath+'global_sigma_'+str(n+1)+'.txt')
    s[:,n] = globalsigma[:,0]
    for i in range(n_genes):
        signal[:,i] = signal[:,i]+output_array[:,i,0] #Right now is just the number of RNAPs transcribing

#Now we need to normalize the signal
#for i in range(n_genes):
#    signal[:,i] =  signal[:,i]/np.max( signal[:,i] )

#And calculate average sigma and std per frame/time
for k in range(frames):
    sigma[:,0] = np.mean( s[k,:] )
    sigma[:,1] = np.std( s[k,:] )

#This t0 is the time at which the system should reach the steady state? or a ~cte global sigma
t0 = 1000 #600#10 mins of equilibration?
#N = len(signal[t0:,0,0])


#----------------------------------------------------------
#TIME TO PLOT
#----------------------------------------------------------
for i in range(n_genes):
    gname = genes_df.iloc[i]['name']
    if gname == 'CDS' or gname == 'tetA':
        continue
    ax.plot( time, signal[:,i], color=colors[i], label=gname )

#ax.plot( time, globalsigma, color='black', label='global')
ax.legend(loc='best')
ax.grid(True)
ax.set_ylabel( 'output')
ax.set_xlabel( 'time (seconds)')
ax.plot( [t0,t0], [-10,10], 'k', lw=10, alpha=0.2 )
ax.set_ylim(0,8)
#ax.set_ylim(0,1.5)

plt.savefig('avg-output.pdf')
