import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from csv import reader
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys

#Description
#--------------------------------------------------------------------------------------------
#I need to plot the sequence and heatmaps of P and G

#Initial conditions
#--------------------------------------------------------------------------------------------
path="../filtering/"
gene_info="../gene.csv"

title = "pbr322_step5"

#Plotting variables
coutput = "profiles_pbr322_step5"
#P_colour = 'Blues_r'
P_colour = 'viridis'
G_colour = 'magma'
my_colour = "blue"
aspectr  = 1#35
width = 7
height = 4

#Drawing genome colours
ah = 115         #height
awidth=40.0      #width of arrow
gene_color='red' #gene color
y_pos = [0, 20, 40, 60, 80, 100]
y_sigma = [-0.1, -0.08, -0.06, -0.04, -0.02, 0.0]

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
w=100 #window size
a=0
for i in range(w-1):   #100 is the size of the window
    sequence_c.append(sequence[a])
    a=a+1

#Let's calculate the content
content = gc_content_window(sequence_c)
content = np.roll(content,int(w/2) ) #It needs to move half of the window so it representes the sourroundins.


#Load and prepare the heatmaps (P and G)
#--------------------------------------------------------------------------------------------
#Info file to read the data
info = np.loadtxt(path+"info.txt")

#I need a test file to allocate the data to the heatmap.
#I assume that there will always be a G_0.txt
test = np.loadtxt( path+"G_0.txt" )

nbp = len(test[:,0])      #number of bp
n_sigma = len(info[:,0])  #number of supercoiling densities
print(nbp, n_sigma)

#Initialize arrays
G_map = np.zeros( ( n_sigma, nbp ) )
P_map = np.zeros( ( n_sigma, nbp ) )
sigma = np.zeros( n_sigma )
bp_ax = test[:,0]

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
fig, axs = plt.subplots(3, figsize=(width,height*3), sharex=True, sharey=False)#, tight_layout=True)
#fig, axs = plt.subplots(3, figsize=(width,height*3), sharex=True, sharey=False, constrained_layout=True)

#Plot sequence
#------------------------------------------------------------------
ax = axs[0]
x = np.arange(1,len(content)+1)
y = content

ax.fill_between(x, y,color=my_colour,alpha=0.5)
ax.set_ylabel(r'GC content (%)')
ax.set_ylim(0,100)
ax.grid(True)
#plt.title(r'$\sigma=$'+c)

#Draw genome----------------------------

#Draw line
line = mlines.Line2D( (0,nbp), (ah, ah), lw=1.5, clip_on=False,color='black' )
ax.add_line(line)

# skip first line i.e. read header first and then iterate over each row od csv as a list
with open(gene_info, 'r') as read_obj:
    csv_reader = reader(read_obj)
    header = next(csv_reader)
    #format is: type, name, start, end, direction
    # Check file as empty
    if header != None:
        # Iterate over each row after the header in the csv
        for row in csv_reader:
            print(row)
            ctype = row[0]
            cname = row[1]
            start = int(row[2]) - 1  #-1 because python counts from 0
            end   = int(row[3]) - 1 
            direc = row[4]

            #Draw genes
            if ctype == "gene":
                obj = mpatches.Arrow( start, ah, end-start, 0, width=awidth, clip_on=False, color=gene_color )
                ax.add_patch(obj)
                #Draw text
                if direc == "+":
                    a = start + abs(end-start)/4
                if direc == "-":
                    a = end + abs(end-start)/4
                ax.text(a, ah+awidth/4, cname, clip_on=False, fontweight='bold')

#Plot G
#------------------------------------------------------------------
ax = axs[1]
im= ax.imshow( G_map, cmap= G_colour, interpolation= 'nearest', aspect=aspectr )

#plt.colorbar(im, shrink=0.7, label=r'$G$ (kcal/mol)', ax=ax, orientation='horizontal')
axins = inset_axes(ax,
                   width="2.5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
fig.colorbar(im, cax=axins, label=r'$G$ (kcal/mol)')# ticks=[1, 2, 3])


#custom y-axis
#y_pos = np.arange( 0, n_sigma, n_sigma/4. )
#ylim = [ y_pos[0], y_pos[-1] ]
#y_sigma = []
#for i in range( len(y_pos) ):
#    j = int( y_pos[i] )
#    y_sigma.append(sigma[j])
ax.set_yticks( y_pos )
ax.set_yticklabels( y_sigma)
ax.set_ylabel(r'$\sigma$')
ax.grid(True)

#Plot P
#------------------------------------------------------------------
ax = axs[2]
im= ax.imshow( P_map, cmap= P_colour, interpolation= 'nearest', aspect=aspectr )
axins = inset_axes(ax,
                   width="2.5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
fig.colorbar(im, cax=axins, label=r'$P$')# ticks=[1, 2, 3])


#plt.colorbar(im, shrink=0.8, label=r'$P$', ax=ax, orientation='horizontal',pad=0.3)
#axins = inset_axes(ax,
#                   width="75%",  # width = 5% of parent_bbox width
#                   height="5%",  # height : 50%
#                   #loc='lower left',
#                   loc='upper center',
#                   bbox_to_anchor=(0.0, 0.2,1, 1),
#                   bbox_transform=ax.transAxes,
#                   borderpad=0,
#                   )
#fig.colorbar(im, cax=axins, orientation='horizontal')# ticks=[1, 2, 3])


ax.set_yticks( y_pos )
ax.set_yticklabels( y_sigma)
ax.set_ylabel(r'$\sigma$')
ax.set_xlabel('bp position')
ax.grid(True)


plt.subplots_adjust(hspace=.05)

plt.savefig(coutput+".pdf")
plt.savefig(coutput+".png")



