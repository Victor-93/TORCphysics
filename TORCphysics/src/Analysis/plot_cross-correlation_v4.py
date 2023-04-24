import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import params

#from scipy import signal


#--------------------------------------------------------------------------
#DESCRIPTION
#--------------------------------------------------------------------------
#I just want to plot the gene system, the transcription signals as blocks,
#and the cross-correlation as done in 2022 Galloway

#params
#---------------------------------------------------------
width = 6
height = 2

circular = True

df = pd.read_csv( 'input.csv', sep='\t')

genome_df = df[ df[ 'type' ].isin(['genome']) ]          #dataframe of genome (to get size)
genome_df = genome_df.reset_index(drop=True)

genes_df  = df[ df[ 'type' ].isin(['gene']) ]            #dataframe of genes (to get orientations...)
genes_df = genes_df.reset_index(drop=True)

n_genes = len(genes_df['start']) #number of genes

o_df  = df[ df[ 'type' ].isin(['NAP']) ]   #dataframe of objects (NAPs )
n_NAPs = len(o_df['start'])           #number of objects bound to the DNA (NAPs)

#Let's load the output
output_array = np.load('results/output_1.npy')
#output_array = np.load('output.npy')

frames = len(output_array[:,0,0] )

time = np.arange(0, frames * params.dt, params.dt)

#Read global sigma
globalsigma = np.loadtxt('results/global_sigma_1.txt')
globalsigma = globalsigma[:,0]


#Better use these ones for the colors of genes...
#Sort them according the input file...
colors = []
colors.append("yellow")
colors.append("black")
colors.append("blue")
colors.append("magenta")
colors.append("red")
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

fig = plt.figure(figsize=(width*2, height*4.5),tight_layout=True)#constrained_layout=True)
gs = fig.add_gridspec(3, 5)
ax1 = fig.add_subplot(gs[:, 0:2], projection='polar')
ax1_l = fig.add_subplot(gs[0, 0:2])
ax2_l = fig.add_subplot(gs[2, 0:2])
ax2 = fig.add_subplot(gs[0, 2:])
ax3 = fig.add_subplot(gs[1, 2:])
ax4 = fig.add_subplot(gs[2, 2:])

#---------------------------------------------------------
#Draw genes
#---------------------------------------------------------
ax = ax1
h=1.5
dh=0.5
hl = 10 #for the legends/labels/info plots
dhl = 1.75#1.5
wl  = .1 #Where genes start to be drawn
dwl  = 1 #genes end at wl+dwl

gene_lw   = 8
object_text = 10 #NAPs, genes, etc...
DNA_lw    = 3
DNA_colour = 'black'

ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

g0 = genome_df.iloc[0]['start']
g1 = genome_df.iloc[0]['end']
nbp = genome_df.iloc[0]['end']
gx = np.array( [g0, g1] )
gy = np.array( [h, h] )
if circular:
    theta = np.arange( 0, 2*np.pi, .1*np.pi)
    radius = 1
    ax.set_rmax( 2)
    gx = radius*np.cos(theta)
    gy = radius*np.sin(theta)
    ax.set_ylim(0,radius*1.2)
else:
    ax.set_xlim(g0-100,g1+100)
    ax.set_ylim(0,2)

for i in range(2):
    if i == 0:
        ax_l = ax1_l
    else:
        ax_l = ax2_l
    ax_l.set_xlim(0,3)
    ax_l.set_ylim(0,hl)
    ax_l.grid(False)
    ax_l.patch.set_alpha(0.0)
    # Hide the right and top spines
    ax_l.spines['right'].set_visible(False)
    ax_l.spines['top'].set_visible(False)
    ax_l.spines['left'].set_visible(False)
    ax_l.spines['bottom'].set_visible(False)
    ax_l.tick_params( left=False, labelleft=False, bottom=False, labelbottom=False, top=False)

xlabels = np.arange(1, nbp+1, int(nbp/8) )
ax.set_xticklabels(xlabels)

#Plot DNA
theta = np.arange( 0, 2*np.pi+.1, 0.01)
r     = np.ones( len(theta) )
ax.plot( theta, r, lw=DNA_lw, color=DNA_colour, zorder=2)

ax.tick_params( labelleft=False, bottom=False,top=False)
#ax.grid(False)

for i in range(n_genes):
    dx = genes_df.iloc[i]['end'] - genes_df.iloc[i]['start']
    x0 = genes_df.iloc[i]['start']
    x1 = genes_df.iloc[i]['end']
    theta0 = 2*np.pi*(x0-1)/nbp
    theta1 = 2*np.pi*(x1-1)/nbp
    if theta0 < theta1:
        theta = np.arange( theta0, theta1, .01  )
    else:
        theta = np.arange( theta0, theta1, -.01  )

    r     = radius*np.ones( len(theta) )

    r     = radius*np.ones( len(theta) )
    
    a0 = np.degrees(theta0)
    a1 = np.degrees(theta1)
    name = genes_df.iloc[i]['name']
    rate = genes_df.iloc[i]['rate']

    ax.plot(theta[:-5], r[:-5], color=colors[i], lw=gene_lw)
    a0 = np.array( [ theta[-2], r[-2]] )
    a1 = np.array( [ theta[-1], r[-1]] )
    arrow = mpatches.FancyArrowPatch( a0, a1, mutation_scale=20, color=colors[i], lw=gene_lw*.75, arrowstyle="simple", zorder=3)
    ax.add_patch(arrow)

    #Labels
    arrow = mpatches.FancyArrowPatch( (wl,hl-dhl*(1+i)), (wl+dwl,hl-dhl*(1+i)), mutation_scale=20, color=colors[i], zorder=3, lw=gene_lw*.75)
    ax1_l.add_patch(arrow)

    #And text
    ax1_l.text(wl+dwl*1.3, hl-dhl*(1+i), name, fontsize=object_text) #name
    ax1_l.text(wl+dwl*2.3, hl-dhl*(1+i), r'$k_0=$'+str(rate), fontsize=object_text) #rate

#Draw naps
for i in range(n_NAPs):
    x0 = o_df.iloc[i]['start']
    name = o_df.iloc[i]['name']
    if circular:
        theta = 2*np.pi*(x0-1)/nbp
        x = radius*np.cos(theta)
        y = radius*np.sin(theta)
        y = radius
        x = theta
        dy = radius*0.15
    else:
        x = x0
        y = h-h/10

    if name == 'ori':
        m = ori_s
        c = ori_c
        ms = ori_ms
    else:
        m = NAP_s
        c = NAP_c
        ms = NAP_ms

    #Plot
    ax.plot( x, y, marker=m, markersize=ms, color = c, mec='k')

    #And labels
    ax2_l.plot( wl+dwl*1.3, hl-dhl*(2.5+i), marker=m, markersize=NAP_ms, color = c, mec='k' )
    ax2_l.text( wl+dwl*2.3,hl-dhl*(2.5+i), name, fontsize=object_text) #name

ax.set_xlabel( 'bp position')

#---------------------------------------------------------
#PLOT SUPERCOILING AND SIGNALS
#----------------------------------------------------------

#Let's First create our signals
signal = np.zeros( ( frames, n_genes ) )

for i in range(n_genes):
    for k in range(frames):
        if output_array[k,i,0] > 0: #If there is one RNAP transcribing
            signal[k,i] = 1

#This t0 is the time at which the cross-correlation takes place.
t0 = 100 #600#10 mins of equilibration?
N = len(signal[t0:])

#This lag is the x axis of the cross correlation plot
lag = np.arange(-N * params.dt * .5, N * params.dt * .5, params.dt)


#Supercoiling
#--------------------------
ax = ax2
for i in range(n_genes):
    if genes_df.iloc[i]['name'] == 'CDS':
        continue
    sigma = output_array[:,i,2]
    ax.plot( time, sigma, color=colors[i] )
ax.plot( time, globalsigma, color='black', label='global')
ax.legend(loc='best')
ax.grid(True)
ax.set_ylabel( 'Supercoiling $\sigma$')
ax.set_xlabel( 'time (seconds)')
ax.plot( [t0,t0], [-10,10], 'k', lw=10, alpha=0.2 )
ax.set_ylim(-.2,.051)


#SIGNAL
#--------------------------
ax = ax3
#Now, let's plot for real
for i in range(n_genes):
    if genes_df.iloc[i]['name'] == 'CDS':
        continue
    ax.fill_between( time, signal[:,i], color=colors[i], alpha=0.5, label=genes_df.iloc[i]['name'] )
ax.grid(True)
ax.set_ylabel( 'Transcription signal')
ax.set_xlabel( 'time (seconds)')

#CROSS CORRELATION
#--------------------------
ax = ax4
acorr = [] #autocorrelations
crossc = [] #cross-correlations with tetA

#Let's first compute auto correlations
for i in range(n_genes):
    a = np.correlate( signal[t0:,i], signal[t0:,i],"same" ) #this is the autocorrelation
    acorr.append(a) #note that it is not normalized

#Now the cross-correlations with tetA gene (in this case is index 0) 
for i in range(n_genes):
    a=np.correlate( signal[t0:,0], signal[t0:,i],"same" )
    crossc.append(a) #It is not normalized

#and PLOT THEM
for i in range(n_genes):
    if genes_df.iloc[i]['name'] == 'CDS':
        continue
    if genes_df.iloc[i]['name'] == 'tetA':
        continue
    if np.max(acorr[i]) > 0:
        a = crossc[i]/np.sqrt( np.max(acorr[0])*np.max(acorr[i]) ) # This normalizes them
    else:
        a = crossc[i] #Because if not, we divide by 0
    print( "max cross-correlation in gene", genes_df.iloc[i]['name'], " at (corr,timelag):", np.max(a), lag[np.argmax(a)] )
    ax.plot( lag, a, color=colors[i] )
#ax.legend(loc='best')
ax.grid(True)
ax.set_xlim(-750,750)
ax.set_ylabel( 'Cross-corr w tetA')
ax.set_xlabel( 'time lag (seconds)')

plt.savefig('cross-correlation.pdf')
plt.savefig('cross-correlation.png')
