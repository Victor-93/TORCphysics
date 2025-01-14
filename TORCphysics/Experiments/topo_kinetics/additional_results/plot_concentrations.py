import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# Plot concentrations
# #Let's try doing it as a heatmap

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'concentrations_plot'

file_in = "concentrations-01-dt"+str(dt)+".pkl"

# We need these concentrations to build matrix
gyrase_concentration_min = 0.0
gyrase_concentration_max = 100.0
gyrase_concentration_step = 2.0
topoI_concentration_min = 0.0
topoI_concentration_max = 100.0
topoI_concentration_step = 2.0

# -----------------------------------
# FIGURE Params
# -----------------------------------
width = 7
height = 3.5
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
experiment_color = 'black'
rec_model_color = 'green'
# poi_model_color = 'green'

# Choose a colormap
#colormap = cm.inferno
colour_maps = ['inferno_r', 'inferno_r', 'Reds', 'Blues']#, 'cool', 'viridis']

gyrase_color = 'blue'
topoI_color = 'red'


# -----------------------------------
# Load input
# -----------------------------------
with open(file_in, 'rb') as file:
    data = pickle.load(file)

# ----------------------------------------------------------------------------------------------------------------------
# Plot - heatmaps
# ----------------------------------------------------------------------------------------------------------------------
column_cases = ['both_Rx', 'both_Sc']  # Both rows should look identical in principle, at least for the cases that reach equilibration
row_cases = ['supercoiling', 'MM_supercoiling', 'topoI', 'gyrase']
# supercoiling and MM_supercoiling should look similar if the model can predict.

# Let's create the array containing concentrations
gyrase_vars = np.arange(gyrase_concentration_min, gyrase_concentration_max + gyrase_concentration_step, gyrase_concentration_step)
topoI_vars = np.arange(topoI_concentration_min, topoI_concentration_max + topoI_concentration_step, topoI_concentration_step)

g_size = len(gyrase_vars)
t_size = len(topoI_vars)

gyrase_pos = np.arange(0, g_size+1, 10)
gyrase_ticks = np.arange(0, gyrase_concentration_max+1, 20)
topoI_pos = np.arange(0, t_size+1, 10)
topoI_ticks = np.arange(0, topoI_concentration_max+1, 20)

# Create figure
fig, axs = plt.subplots(len(row_cases), len(column_cases), figsize=(width*len(column_cases), len(row_cases) * height), tight_layout=True)

relative_index_list = [] # I will add cases for which we have both enzymes present

for i, col_case in enumerate(column_cases):
    for j, row_case in enumerate(row_cases):

        ax = axs[j, i]

        title = row_case + ' for ' + col_case
        ax.set_title(title, fontsize=title_size)

        # Building matrix
        cmatrix = np.zeros([g_size, t_size])

        # Go through data to extract matrix values
        s = 0
        for gyrase_i, gyrase_concentration in enumerate(gyrase_vars):
            for topoI_i, topoI_concentration in enumerate(topoI_vars):

                if gyrase_concentration >0 and topoI_concentration >0:
                    ratio = gyrase_concentration/topoI_concentration
                else:
                    ratio = 0.0

#                print(gyrase_i, topoI_i, s, col_case, row_case )
                if row_case == 'MM_supercoiling':
                    a = data[s][col_case][row_case][-1]
                else:
                    a = data[s][col_case][row_case][0][-1]

                if 'supercoiling' in row_case:
                    if a >= -0.08 and a <= -0.04 and ratio > 0.0  and ratio<=5 and gyrase_concentration <= 80.0:
                        a=a
                        relative_index_list.append(s)
                    #                        if ratio > 0 and ratio <= 5:
#                            relative_index_list.append(s)

#                        if 'supercoiling' == row_case:
#                            if (gyrase_concentration >= 20.0 and gyrase_concentration <=80.0 and
#                                    topoI_concentration >= 10.0 and topoI_concentration <=40.0):
#                            if ratio > 0 and ratio <=5:
#                                relative_index_list.append(s)
                    else:
                        a = 0.0
 #               print(a)
                cmatrix[gyrase_i, topoI_i]  = a
                s=s+1

        im  = ax.imshow( cmatrix, cmap= colour_maps[j], interpolation= 'nearest')#, vmin=-0.11, vmax=-0.02)#, aspect=aspectr)

        axins = inset_axes(ax,
                       width="2.5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
        fig.colorbar(im, cax=axins)#, label=labels[i])# ticks=[1, 2, 3])

 #       plot_genes(ax, y_genes)

        ax.set_yticks( gyrase_pos )
        ax.set_yticklabels( gyrase_ticks)
        ax.set_xticks(topoI_pos)
        ax.set_xticklabels(topoI_ticks)
        ax.set_ylabel(r'Gyrase concentration (nM)')
        ax.set_xlabel(r'TopoI concentration (nM)')
#       ax.set_xlabel('bp position')
 #       ax.grid(True)

plt.show()

# ----------------------------------------------------------------------------------------------------------------------
# Plot - relatives concetrations
# ----------------------------------------------------------------------------------------------------------------------
row_cases = ['supercoiling', 'MM_supercoiling', 'topoI', 'gyrase']

# Create figure
fig, axs = plt.subplots(len(row_cases), len(column_cases), figsize=(width*len(column_cases), len(row_cases) * height), tight_layout=True)

for i, col_case in enumerate(column_cases):
    for j, row_case in enumerate(row_cases):

        ax = axs[j, i]

        title = row_case + ' for ' + col_case
        ax.set_title(title, fontsize=title_size)

        x = []
        y = []
        ys = []

        for s in relative_index_list:
            rel = data[s]['gyrase_concentration']/data[s]['topoI_concentration']

            x.append(rel)

            if row_case == 'MM_supercoiling':
                y.append( data[s][col_case][row_case][-1])
            else:
                y.append( data[s][col_case][row_case][0][-1])
                ys.append( data[s][col_case][row_case][1][-1])

        ax.plot(x, y, 'o', color='red')

        ax.set_ylabel(row_case, fontsize=xlabel_size)
        ax.set_xlabel(r'Gyrase / Topo I ratio', fontsize=xlabel_size)
        ax.grid(True)

plt.show()
