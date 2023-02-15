import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import pandas as pd
import params
import statistical_model as sm
import sys

#params
#---------------------------------------------------------
width = 6
height = 4

x = np.arange(-.121,.04,.001)
df = pd.read_csv( 'input.csv', sep='\t')

genes_df  = df[ df[ 'type' ].isin(['gene']) ]            #dataframe of genes (to get orientations...)
genes_df = genes_df.reset_index(drop=True)

n_genes = len(genes_df['start']) #number of genes


colors = []
colors.append("green")
colors.append("brown")
colors.append("blue")
colors.append("magenta")
colors.append("yellow")
colors.append("purple")
colors.append("cyan")
colors.append("black")

#plot
#---------------------------------------------------------

fig, ax = plt.subplots(1, figsize=(width, height), tight_layout=True)

for i in range(n_genes):
    if genes_df.iloc[i]['name'] == 'tetA':
        rate = np.zeros( len(x) )
        rate[:] = genes_df.iloc[i]['rate']
    else:
        rate = sm.promoter_curve( genes_df.iloc[i]['rate'], x )
    ax.plot(x,rate, color=colors[i], label=genes_df.iloc[i]['name'] )

ax.legend(loc='best')
ax.grid(True)
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'Initiation rate ($s^{-1}$)')
ax.set_title('promoter response curve')

plt.savefig("promoter_curves.png")
