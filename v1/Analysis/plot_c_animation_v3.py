import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import pandas as pd
import params
from datetime import datetime
import statistical_model as sm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patheffects as pe
import sys


#--------------------------------------------------------------------------
#DESCRIPTION
#--------------------------------------------------------------------------
#I just want to create a cute animation FOR CIRCULAR DNA

#INPUTS
#--------------------------------------------------------------------------

coutput="circular_animation"

frames = 3000
width = 12
height = 6

#files
posfile='results/positions_1.txt'
sigmafile='results/supercoiling_1.txt'
namefile='results/object_1.txt'
cinput = "input.csv"

globalsigma = np.loadtxt('results/global_sigma_1.txt')

sfile = open(sigmafile, "r")
pfile = open(posfile, "r")
nfile = open(namefile, "r")

x = np.arange(-.15,.15,.001)

topo = sm.topoisomerase_activity(x)
gyra = sm.gyrase_activity(x)

#---------------------------------------------------------------------------------------------------------------------
#PARAMETERS
#---------------------------------------------------------------------------------------------------------------------
#All parameters are already in the params module, but I prefer to have them here with more simple names:
v0     = params.v0
w0     = params.w0
gamma  = params.gamma
dt     = params.dt
ds     = 1 #This is the interval of frames we are going to skip to speed up the animation process

#------------------------------------------------------------
#FUNCTIONS
#------------------------------------------------------------

#Reads a line
#------------------------------------------------------------
def build_RNAP_arrays():
    #We need to get:
    #The position of the domains x[j]
    #their height                y[j]
    #The supercoiling value  sigma[j]
    #And the position of RNAPs   rx[s], ry[s]
    pos  = pfile.readline().split()
    sup  = sfile.readline().split()
    name = nfile.readline().split()

    n = len(pos)
    x = pos
    y = h*np.ones(n)
    x = np.zeros(n)
    sigma = np.zeros(n)
    rx = -2000*np.ones(n)
    ry = y#-100*np.ones(n)
    s = 0
    for j in range(n):
        x[j] = pos[j]
        sigma[j] = sup[j]
        if name[j] == 'RNAP':
            rx[j] = pos[j]

    return x,y,rx,ry,sigma

#Animates the curves
#------------------------------------------------------------
def animate(i):
    RNAP_plot.set_data( RNAP_theta[i], RNAP_radius[i])
    
    #Let's take the values we need...
    n = n_ob[i]
    theta = stg_theta[i]
    rar   = sigma_radius[i]
    c     = sigma_color[i]
    t     = topo_radius[i]
    g     = gyra_radius[i]
    
    for j in range(10):
        lines[j].set_data( [1,1.02], [100,100]  )
        lines_topo[j].set_data( [1,1.02], [100,100]  )
        lines_gyra[j].set_data( [1,1.02], [100,100]  )

        lines[j].set_linewidth(.2)
        lines_topo[j].set_linewidth(.2)
        lines_gyra[j].set_linewidth(.2)

    for j in range(n): #And plot them as we per domain j

        lines[j].set_data( theta[j], rar[j])
        lines[j].set_color(c[j])
        lines_topo[j].set_data( theta[j], t[j])
        lines_gyra[j].set_data( theta[j], g[j])

        lines[j].set_linewidth(sigma_lw)
        lines_topo[j].set_linewidth(topo_lw)
        lines_gyra[j].set_linewidth(topo_lw)


    #for the time
    time_text.set_text( my_time[i] )

    
    #return lines, RNAP_plot, time_text#, lines_topo, lines_gyra#, lines_tg
    return RNAP_plot, lines_topo, lines_gyra, time_text


#PLOTTING VARIABLES
#--------------------------------------------------------------------------
#colours
gene_colour = 'green'
DNA_colour = 'black'
sigma_colour = 'red'
ori_c  = "yellow"
NAP_c  = "gray"
RNAP_c  = "white"

#Better use these ones for the colors of genes...
#Sort them according the input file...
colors = []
#colors.append("purple")
colors.append("yellow")
colors.append("black")
#colors.append("pink")
colors.append("blue")
colors.append("magenta")
colors.append("red")
colors.append("green")
colors.append("cyan")
colors.append("black")

#Sizes
NAP_ms = 8
ori_ms = 15
ori_s  = 'p'
NAP_s  = "s"
RNAP_ms = 20
RNAP_s  = "o"
DNA_lw    = 5#6#12
gene_lw   = 10#5
sigma_lw  = 20
topo_lw   = 4

#text size
slabel = 15
object_text = 10 #NAPs, genes, etc...

#Color maps
cmap = cm.get_cmap('RdBu') #This one is for supercoiling
#cmap = cm.get_cmap('seismic_r') #This one is for supercoiling
#cmap = cm.get_cmap('cividis') #This one is for supercoiling
#cmap = cm.get_cmap('RdBu_r') #This one is for supercoiling
#cmap = cm.get_cmap('cool') #This one is for supercoiling
sigma_min = -.06#np.min( globalsigma[:,0] )
sigma_max = .06#.02#np.max( globalsigma[:,0] )
print(sigma_min, sigma_max)
#sigma_min = -.8
#sigma_max = .2
sigma_m   = sigma_max - sigma_min #This defines the width...

#-----------------------------------
#Read genome
#-----------------------------------

#Sizes of plot
#-----------------------------------
topo_fact= 1#60/dt #Factor by which topo is going to be scaled
topo_m = np.max(topo)*topo_fact #Topo adds positive supercoils when sigma is negative
gyra_m = np.min(gyra)*topo_fact #And gyra adds negative supercoils when sigma is positive
topo_w = topo_m - gyra_m
topo_c = "black"#cmap(1.0)
gyra_c = "gray"#cmap(0.0)


ylima = 0
if topo_m > abs(gyra_m):
    ylimb = topo_m + topo_w*.3 
else:
    ylimb = abs(gyra_m) + topo_w*.3 

h=ylimb-topo_w*.1#topo_m + topo_w*.3#175#2#1.5
dh=h/4#0.5
hl = 10 #for the legends/labels/info plots
dhl = 1
wl  = 1.2 #Where genes start to be drawn
dwl  = 1 #genes end at wl+dwl

#Read genome
#-----------------------------------
df = pd.read_csv(cinput,sep='\t')
genome_df = df[ df[ 'type' ].isin(['genome']) ]          #dataframe of genome (to get size)
genome_df = genome_df.reset_index(drop=True)
genes_df  = df[ df[ 'type' ].isin(['gene']) ]            #dataframe of genes (to get orientations...)
genes_df = genes_df.reset_index(drop=True)
n_genes = len(genes_df['start']) #number of genes

o_df  = df[ df[ 'type' ].isin(['NAP']) ]   #dataframe of objects (NAPs )
n_NAPs = len(o_df['start'])           #number of objects bound to the DNA (NAPs)

#More sizes
#-----------------------------------
g0 = genome_df.iloc[0]['start']
g1 = genome_df.iloc[0]['end']
gx = np.array( [g0, g1] )
gy = np.array( [h, h] )
nbp = genome_df.iloc[0]['end']

#--------------------------------------------------------------------------
#PLOT
#--------------------------------------------------------------------------
fig = plt.figure(figsize=(width, height),tight_layout=True)#constrained_layout=True)
gs = fig.add_gridspec(1, 6)
ax = fig.add_subplot(gs[0, 0:4], projection='polar')  #this one is for plotting
ax_l = fig.add_subplot(gs[0, 4:]) # This one is for legends

#Sizes and limits
#-----------------------------------
ax.set_theta_zero_location("N") #This aligns the 0 location at the north position
ax.set_theta_direction(-1)      #The orientation of the gene

theta = np.arange( 0, 2*np.pi, .1*np.pi)
radius = h#1.2
gx = radius*np.cos(theta)
gy = radius*np.sin(theta)
ax.set_ylim(ylima,ylimb)
y = np.arange(0, ylimb-topo_w*.3 , topo_w/5)
#ylabel = str( y*60/dt)
ax.set_yticks( y )
ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax_l.set_xlim(0,3)
ax_l.set_ylim(0,hl)

xlabels = np.arange(1, nbp+1, int(nbp/8) )
ax.set_xticklabels(xlabels)

#labels and all that
#-----------------------------------
ax.grid(True, zorder=1)
ax_l.grid(False)

# Hide the right and top spines
ax_l.spines['right'].set_visible(False)
ax_l.spines['top'].set_visible(False)
ax_l.spines['left'].set_visible(False)
ax_l.spines['bottom'].set_visible(False)

ax_l.tick_params( left=False, labelleft=False, bottom=False, labelbottom=False, top=False)

#-----------------------------------
#Now draw genes
#-----------------------------------
n_genes = len(genes_df['start'])
theta = np.arange( 0, 2*np.pi+.1, 0.01)
r     = radius*np.ones( len(theta) )
ax.plot( theta, r, lw=DNA_lw, color=DNA_colour, zorder=2)

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

    a0 = np.degrees(theta0)
    a1 = np.degrees(theta1)
    name = genes_df.iloc[i]['name']
    rate = genes_df.iloc[i]['rate']

    ax.plot(theta[:-5], r[:-5], color=colors[i], lw=gene_lw, path_effects=[pe.Stroke(linewidth=gene_lw+3, foreground='k'), pe.Normal()])
    a0 = np.array( [ theta[-2], r[-2]] )
    a1 = np.array( [ theta[-1], r[-1]] )
    arrow = mpatches.FancyArrowPatch( a0, a1, mutation_scale=20, color=colors[i], lw=gene_lw*.75, arrowstyle="simple", zorder=3, path_effects=[pe.Stroke(linewidth=gene_lw, foreground='k'), pe.Normal()])
    ax.add_patch(arrow)

    #Add labels
    arrow = mpatches.FancyArrowPatch( (wl,hl-dhl*(1+i)), (wl+dwl,hl-dhl*(1+i)), mutation_scale=25, color=colors[i], zorder=3, lw=gene_lw)
    ax_l.add_patch(arrow)

    #And the text
    ax_l.text(wl+dwl*1.3, hl-dhl*(1+i), name, fontsize=object_text) #name
    ax_l.text(wl+dwl*1.3, hl-dhl*(1.3+i), r'$k_0=$'+str(rate), fontsize=object_text) #rate
    
#-----------------------------------
#Draw NAPs
#-----------------------------------
for i in range(n_NAPs):
    x0 = o_df.iloc[i]['start']
    name = o_df.iloc[i]['name']
    theta = 2*np.pi*(x0-1)/nbp
    y = radius
    x = theta
    dy = radius*0.15

    if name == 'ori':
        m = ori_s
        c = ori_c
        ms = ori_ms
    else:
        m = NAP_s
        c = NAP_c
        ms = NAP_ms

    #Plot
    ax.plot( x, y, marker=m, markersize=ms, color = c, mec='k', zorder=4 )

    #And Labels
    ax_l.plot( wl+dwl*.5, hl-dhl*(1+i+n_genes), marker=m, markersize=NAP_ms, color = c, mec='k' )
    ax_l.text( wl+dwl*1.3,hl-dhl*(1.1+i+n_genes), name, fontsize=object_text) #name

#Add RNAP label
ax_l.plot( wl+dwl*.5, hl-dhl*(1+n_genes+n_NAPs), marker=RNAP_s, markersize=RNAP_ms, color = RNAP_c, mec='k' )
ax_l.text( wl+dwl*1.3, hl-dhl*(1.1+n_genes+n_NAPs), "RNAP", fontsize=object_text) #name

#-----------------------------------
#Draw color bar of sigma
#-----------------------------------
rw = 1
rh = .1
y0 = 0.2
x0 = wl-rw*2
x = x0
for i in range(101):   #DRAW BAR
    y = y0 + i*rh
    a = i/100
    c = cmap( a )

    rec = mpatches.Rectangle( (x,y), rw, rh, color = c )
    ax_l.add_patch( rec)

#DRAW TICKS
y = np.arange( y0, y0+100*rh+1, 100*rh/5 )
sup = np.arange( sigma_min, sigma_max+sigma_m, sigma_m/5 )
yn = len(y)
x0 = wl-rw
for i in range(yn):

    z = sup[i]
    ax_l.text( x0, y[i]-y0, f'{z:.3f}')

#And write label
ax_l.text( wl-rw*.4, hl/2.25, r'local $\sigma$', fontsize=slabel, rotation = 90  )

#-----------------------------------
#Draw topoisomerase legends
#-----------------------------------
#Topo I
x = [ wl, wl+dwl ]
y = [ hl-dhl*(2+n_genes+n_NAPs), hl-dhl*(2+n_genes+n_NAPs) ]
ax_l.plot(  x, y, "-", c = topo_c, lw=topo_lw) 
ax_l.text( wl+dwl*1.3, hl-dhl*(2+n_genes+n_NAPs)+.1, "Topo I", fontsize=object_text)

#Gyrase
y = [ hl-dhl*(3+n_genes+n_NAPs), hl-dhl*(3+n_genes+n_NAPs) ]
ax_l.plot(  x, y, "-", c = gyra_c, lw=topo_lw) 
ax_l.text( wl+dwl*1.3, hl-dhl*(3+n_genes+n_NAPs)+.1, "Gyrase", fontsize=object_text)

#-----------------------------------
#THE ANIMATION
#-----------------------------------

#BEFORE ANIMATING, LET'S BUILD THE ARRAYS 
#-----------------------------------
RNAP_theta = []
RNAP_radius = []
stg_theta = []   #Theta for the domains
sigma_radius = []
sigma_color  = []
topo_radius = []
gyra_radius = []
n_ob = []
my_time = []

o = 2*np.pi*RNAP_ms*2.5/nbp #We need this to adjust the domains in the middle of the RNAPs
s = 0 #Let's reduce the number of frames
for k in range(frames):
    if s == ds:
        s = 0
    s = s+1
    x,y,rx,ry,sigma = build_RNAP_arrays()

    if s != 1 and s < ds+1:
        continue
    n = len(x)
    n_ob.append(n)
    a = []
    b = []
    for j in range(n):
        theta = 2*np.pi*( rx[j] - 1 )/nbp
        a.append(theta)
        if rx[j] <= -1000:
            b.append( 1000)
        else:
            b.append(radius)

    #DATA FOR RNAP
    RNAP_theta.append(a)
    RNAP_radius.append(b)

    theta = []
    rar   = []
    to    = []
    gy    = []
    c     = []

    for j in range(n):
        t = sm.topoisomerase_activity( sigma[j])*topo_fact
        g = abs( sm.gyrase_activity( sigma[j]) )*topo_fact #Abs because it goes in the range [0,a] 

        if j < n-1:
            th0 = 2*np.pi* ( x[j] - 1)/nbp
            th1 = 2*np.pi* ( x[j+1] -1 )/nbp
        if j == n-1:
            th0 = 2*np.pi* ( x[j] - 1)/nbp
            th1= 2*np.pi* ( g1 - 1)/nbp

        th = np.arange( th0, th1, 0.1) 
        r = radius*np.ones(len(th))
        t = t*np.ones(len(th))
        g = g*np.ones(len(th))
        theta.append(th+o)
        rar.append(r)
        to.append(t)
        gy.append(g)
        a = (sigma[j]-sigma_min)/sigma_m
        c.append( cmap( a ) )

    stg_theta.append(theta)
    sigma_radius.append(rar)
    sigma_color.append( c )
    topo_radius.append(to)
    gyra_radius.append(gy)

    #And the time
    ttt = datetime.fromtimestamp( k*dt-3600).strftime('%H:%M:%S')
    my_time.append( 'time='+ttt)

new_frames = len(my_time)

print( "Writing ", new_frames, "frames")
#sys.exit()

#INITIALIZE
#-----------------------------------
RNAP_plot, = ax.plot( [], [], ms=RNAP_ms, c=RNAP_c, marker=RNAP_s, zorder=4, mec='k', lw=0)

lines = [ax.plot( [], [], c=sigma_colour, lw=sigma_lw, zorder=1)[0] for _ in range(20)] #This plots supercoiling

lines_topo = [ax.plot( [], [], "-",  c=topo_c, lw=topo_lw, alpha=1.0)[0] for _ in range(20)] #This is topo I activity
lines_gyra = [ax.plot( [], [], "-", c=gyra_c, lw=topo_lw, alpha=1.0)[0] for _ in range(20)] #This is gyrase activity

#Prepare time_text which will be used for writing the variable time
time_text = ax.text(0.0, 0.98, '', transform=ax.transAxes)

#ANIMATE
#-----------------------------------
ani = animation.FuncAnimation(
    fig, animate, interval=600, frames=500)#new_frames)#, blit=True, frames=200, repeat=True)
    #fig, animate, interval=400, frames=50, save_count=0, cache_frame_data=0)#int(frames/10))
    #fig, animate, frames=50)#int(frames/10))

#SAVE OR SHOW
#-----------------------------------
#plt.show()

writervideo = animation.FFMpegWriter(fps=5, bitrate=400)#1800)
ani.save(coutput+".mp4", writer=writervideo)
#ani.save('myAnimation.gif', writer='imagemagick', fps=30)

