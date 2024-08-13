import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import matplotlib.animation as animation
import numpy as np
import matplotlib.patches as mpatches
from TORCphysics import analysis as an
from TORCphysics import effect_model as em
from datetime import datetime
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import os
from TORCphysics.src import enzyme_shapes_dir
import random

# This module is produced to make one's life easier and be able to quickly process results.
# With this module you can create animations,

# Figure parmameters
width = 10
height = 6

# colors
enzyme_colors = {'RNAP': 'white', 'IHF': 'yellow', 'FIS': 'red', 'lacI': 'green', 'lacI_bridge': 'red', 'ori': 'silver',
                 'topoI': 'red',
                 'gyrase': 'cyan'}

gene_colour = '#4a86e8ff'  # 'blue'
DNA_colour = 'black'  # '#999999ff' #'gray'
sigma_colour = 'red'
torque_colour = 'purple'
topoI_colour = 'red'
gyrase_colour = 'cyan'

# Sizes
enzyme_sizes = {'RNAP': 300, 'IHF': 500, 'FIS': 500, 'lacI': 250, 'lacI_bridge': 250, 'ori': 500, 'topoI': 500,
                'gyrase': 500}
DNA_lw = 12
gene_lw = 5
sigma_lw = 5

# Shapes
enzyme_shapes = {'RNAP': 'o', 'IHF': 'o', 'FIS': 'o', 'lacI': 's', 'lacI_bridge': 's', 'ori': 'h', 'topoI': 'X',
                 'gyrase': 'X'}

# text size
slabel = 20
object_text = 15  # NAPs, genes, etc...

# Ranges
sigma_a = -.3
sigma_b = .2


# Path to molecule/enzymes/effectors png files
#effector_path = '../enzyme_shapes/'


# TODO:
#  1.- Plots
#  1.1- Promoter response curve - DONE
#  1.2- Topoisomerase activity curves - DONE
#  1.3- Signal profiles - all or by type (optional) - can only be rectangular pulses?
#  1.4- Supercoiling profiles
#  1.5- Cross-correlation  - should find maximums and time-lags. print them. (optional)
#  1.6- Steady-state - should extract rates and plot them
#  1.7- Environment plots - mRNA counts
#  2.- Representations
#  2.1- Linear
#  2.2- Circular
#  3.- Animations
#  3.1- Linear
#  3.2- Circular
# TODO: Maybe one wants to include specific sites names? instead of ignoring?
# TODO: You need to test these functions
# TODO: Check how these functions change with the stochastic topo binding

def ax_params(axis, xl, yl, grid, legend):
    axis.grid(grid)
    axis.set_ylabel(yl)
    axis.set_xlabel(xl)
    if legend:
        axis.legend(loc='best')


# PLOTTING CALCULATIONS FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
# The following functions plot certain calculations using the analysis module.
# These functions are provided as a tool to quickly analyse data, so the circuit is an input and it analyses it
# according its sites, enzymes and environment.
# Some functions may have the following inputs:
# axs = axes to plot
# colors = is a dictionary of colors in the form {site_name:color}
# site_type = If you want to calculate/plot by certain site_type, e.g. 'gene'
# labels = If you want the function to plot the x,y labels & grid
# **kwargs = Additional plot arguments. Be careful when indicating the colors in the **kwargs when using colors=dict
# [fa,fb] = frame ranges in which certain quantities are measured, e.g. cross-correlations, steady state curve...
# ignore = List of site names to ignore


# Plots the steady state curve: log( sum of initiation events / time ) vs time
def plot_steady_state_initiation_curve(my_circuit, sites_df, axs=None, ignore=None,
                                       fa=0, fb=-1, colors=None, site_type=None, labels=True, **kwargs):
    if axs is None:
        axs = plt.gca()
    time = np.arange(0, my_circuit.dt * (my_circuit.frames + 1), my_circuit.dt)
    if site_type is None:
        curves, rates, names = an.initiation_rates(sites_df, time, fa, fb)
    else:
        curves, rates, names = an.initiation_rates_by_type(sites_df, site_type, time, fa, fb)

    for i, name in enumerate(names):
        curve = curves[i]
        rate = rates[i]
        if ignore is not None:
            if name in ignore:
                continue
        my_label = name + f' k={rate:.4f}s'
        if colors is not None:
            axs.plot(time, curve, color=colors[name], label=my_label, alpha=0.5, **kwargs)
        else:
            axs.plot(time, curve, label=my_label, alpha=0.5, **kwargs)
    if labels:
        ax_params(axis=axs, xl='time (seconds)', yl=r'$\log$($\sum$initiation)/time', grid=True, legend=True)


# TODO: Cross-correlation with no specific site? That might be a mess
# Plot cross-correlations between sites with respect another site with name ref_name
# ref_name is the name of the site that you want to calculate the cross-correlation with. It has to have the same type
# as site_type in case you specify.
def plot_cross_correlation_with_site(my_circuit, sites_df, ref_name, axs=None, ignore=None,
                                     fa=0, fb=-1, colors=None, site_type=None, labels=True, **kwargs):
    if axs is None:
        axs = plt.gca()
    if site_type is None:
        signals, names = an.build_signals(sites_df)
    else:
        signals, names = an.build_signal_by_type(sites_df, site_type)
    signals_t0 = []
    for signal in signals:
        signals_t0.append(signal[fa:fb])
    for i, name in enumerate(names):
        if name == ref_name:
            index = i
    cross, lag = an.cross_correlation_hmatrix(signals_t0, my_circuit.dt)
    maxlag = []
    j = -1
    for i, name in enumerate(names):
        if name == ref_name:
            continue
        if ignore is not None:
            if name in ignore:
                continue
        j += 1
        # We need to find the maximum correlation write it
        maxlag.append(lag[np.argmax(cross[index, i, :])])
        my_label = name + f' lag={maxlag[j]:.2f}s'
        if colors is not None:
            axs.plot(lag, cross[index, i, :], color=colors[name], label=my_label, **kwargs)
        else:
            axs.plot(lag, cross[index, i, :], label=my_label, **kwargs)
    if labels:
        ax_params(axis=axs, xl='time lag (seconds)', yl='Cross-correlation w ' + ref_name, grid=True, legend=True)
    axs.set_xlim(-200, 200)


# Plots supercoiling profiles at the sites and global
def plot_supercoiling_profiles(my_circuit, sites_df, axs=None, ignore=None, colors=None, site_type=None,
                               only_global=False, labels=True, **kwargs):
    if axs is None:
        axs = plt.gca()
    time = np.arange(0, my_circuit.dt * (my_circuit.frames + 1), my_circuit.dt)
    # Let's plot the global superhelical density
    mask = sites_df['type'] == 'circuit'  # This one contains global superhelical density
    global_sigma = sites_df[mask]['superhelical'].to_numpy()
    axs.plot(time, global_sigma, color='black', label='global')  # and the global
    if not only_global:
        # get names
        if site_type is None:
            mask = sites_df['type'] != 'circuit'
            # Let's filter the sites names
            names = sites_df[mask].drop_duplicates(subset='name')['name']
        else:
            mask = sites_df['type'] == site_type
            names = sites_df[mask].drop_duplicates(subset='name')['name']
        # And plot the superhelical density at sites
        for i, name in enumerate(names):
            if ignore is not None:
                if name in ignore:
                    continue
            mask = sites_df['name'] == name
            superhelical = sites_df[mask]['superhelical'].to_numpy()
            if colors is not None:
                axs.plot(time, superhelical, color=colors[name], label=name, alpha=0.5, **kwargs)
            else:
                axs.plot(time, superhelical, label=name, alpha=0.5, **kwargs)

    if labels:
        ax_params(axis=axs, xl='time (seconds)', yl='Supercoiling at site', grid=True, legend=True)


# This one plots the signal profiles.
def plot_signal_profiles(my_circuit, sites_df, axs=None, ignore=None, colors=None, site_type=None,
                         labels=True, **kwargs):
    if axs is None:
        axs = plt.gca()
    if site_type is None:
        signals, names = an.build_signals(sites_df)
    else:
        signals, names = an.build_signal_by_type(sites_df, site_type)
    time = np.arange(0, my_circuit.dt * len(signals[0]), my_circuit.dt)
    y_0s = time * 0.0
    for i, signal in enumerate(signals):
        name = names[i]
        if ignore is not None:
            if name in ignore:
                continue
        if colors is not None:
            axs.plot(time, signal, color=colors[name], label=names[i], alpha=0.5, **kwargs)
            axs.fill_between(time, signal, y_0s, color=colors[name], alpha=0.5, **kwargs)
        else:
            axs.plot(time, signal, label=names[i], alpha=0.5, **kwargs)
    if labels:
        ax_params(axis=axs, xl='time (seconds)', yl='Transcription signal', grid=True, legend=True)


# Plots the site responses (rates) modulated by supercoiling
def plot_site_response_curves(my_circuit, axs=None, ignore=None, colors=None, site_type=None, labels=True, **kwargs):
    if axs is None:
        axs = plt.gca()
    i = -1
    for site in my_circuit.site_list:
        if ignore is not None:
            if site.name in ignore:
                continue
        if site_type is None or site.site_type == site_type:
            i += 1
            environment = [environment for environment in my_circuit.environmental_list
                           if environment.site_type == site.site_type]
            if not environment:
                continue
            else:
                environment = environment[0]
            if site.name.isdigit() and 'DNA' in environment.site_type:  # skip non-specific binding proteins
                continue

            rate, x = an.site_activity_curves(site, environment)
            if colors is not None:
                axs.plot(x, rate, color=colors[site.name], label=site.name, **kwargs)
            else:
                axs.plot(x, rate, label=site.name, **kwargs)

    if labels:
        ax_params(axis=axs, xl=r'$\sigma$', yl=r'Initiation rate ($s^{-1}$)', grid=True, legend=True)


# Plots the topoisomerase activity curves of a continuum model
def plot_topoisomerase_activity_curves_continuum(my_circuit, axs=None, ignore=None, labels=True, **kwargs):
    if axs is None:
        axs = plt.gca()
    i = -1
    for environmental in my_circuit.environmental_list:
        if environmental.enzyme_type == 'topo':
            i += 1
            topo_curve, x = an.topoisomerase_activity_curves_continuum(environmental, dt=my_circuit.dt)
            axs.plot(x, topo_curve, label=environmental.name, **kwargs)
            if i == 0:
                topo_sum = np.zeros_like(topo_curve)
            topo_sum += topo_curve
    axs.plot(x, topo_sum, color='black', label='sum', **kwargs)
    if labels:
        ax_params(axis=axs, xl=r'$\sigma$', yl=r'$\sigma$ removed per timestep', grid=True, legend=True)


# PLOTTING REPRESENTATIONS FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------

# PLOTTING ANIMATIONS FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------

def create_animation_linear(my_circuit, sites_df, enzymes_df, output, out_format,
                            site_type=None, site_colours=None, plot_torques=False):
    # plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'

    output_file = output + out_format
    h = 1.5
    dh = 0.5
    gx = np.array([0, my_circuit.size])
    gy = np.array([h, h])

    torque_a = -200
    torque_b = 200

    if plot_torques:
        fig, ax = plt.subplots(3, figsize=(width, height * 1.5),
                               gridspec_kw={'height_ratios': [1, 1.5, 1.5], 'hspace': 0.325})
    else:
        fig, ax = plt.subplots(2, figsize=(width, height), gridspec_kw={'height_ratios': [1, 2], 'hspace': 0.2})

    # Sizes
    # -----------------------------------
    ax[0].set_xlim(- 100, my_circuit.size + 100)
    ax[0].set_ylim(0, 2)
    ax[1].set_xlim(- 100, my_circuit.size + 100)
    ax[1].set_ylim(sigma_a, sigma_b)
    if plot_torques:
        ax[2].set_xlim(- 100, my_circuit.size + 100)
        ax[2].set_ylim(torque_a, torque_b)

    # labels and all that
    # -----------------------------------
    ax[0].grid(axis='x', zorder=1)
    ax[0].tick_params(labelleft=False, bottom=False, top=False)
    ax[1].grid(True, zorder=1)
    ax[1].set_xlabel("DNA (bp)", fontsize=slabel)
    ax[1].set_ylabel(r"$\sigma$", fontsize=slabel)

    if plot_torques:
        ax[2].grid(True, zorder=1)
        ax[2].set_xlabel("DNA (bp)", fontsize=slabel)
        ax[2].set_ylabel(r"Torque (pN nm)", fontsize=slabel)

    # -----------------------------------
    # draw DNA
    # -----------------------------------
    ax[0].plot(gx, gy, lw=DNA_lw, color=DNA_colour, zorder=2)
    # -----------------------------------
    # Now draw genes
    # -----------------------------------
    n_sites = len(my_circuit.site_list)  # - 2
    for i, site in enumerate(my_circuit.site_list):
        if site.site_type == 'EXT':
            continue
        if site_type is not None and site.site_type != site_type:  # Only site types
            continue
        if 'DNA_' in site.site_type:  # I'm trying to avoid bare DNA binding sites
            continue
        # print(site.name)
        # sitess = [site for site in my_circuit.site_list if site.site_type == 'gene']
        x1 = site.end
        x0 = site.start
        dx = x1 - x0
        name = site.name
        if site_colours is not None:
            arrow = mpatches.FancyArrowPatch((x0, h), (x1, h), mutation_scale=25, color=site_colours[name], zorder=3,
                                             lw=gene_lw)
        else:
            arrow = mpatches.FancyArrowPatch((x0, h), (x1, h), mutation_scale=25, color=gene_colour, zorder=3,
                                             lw=gene_lw)
        ax[0].add_patch(arrow)
        if x0 < x1:
            a = x0 + abs(dx / 2)
        else:
            a = x1 + abs(dx / 2)
        ax[0].text(a, h - dh, name, fontsize=object_text)

    # -----------------------------------
    # THE ANIMATION
    # -----------------------------------

    # Prepare data
    # -----------------------------------
    xl = []  # position in x
    yl = []  # in y
    sl = []  # size
    cl = []  # colour
    ml = []  # marker
    sigma = []  # superhelical
    torque = []  # torque
    my_time = []
    mask = sites_df['type'] == 'circuit'
    n_enzymes_df = sites_df[mask]

    l = -1
    for k in range(my_circuit.frames):
        n_enz = n_enzymes_df.iloc[k]['#enzymes'] + 2  # 2 for the EXT
        x = []
        y = []
        s = []
        c = []
        m = []
        sig = []
        tor = []
        for i in range(n_enz):
            l = l + 1
            name = enzymes_df.iloc[l]['name']
            x.append(enzymes_df.iloc[l]['position'])
            y.append(h)
            sig.append(enzymes_df.iloc[l]['superhelical'])
            tor.append(em.Marko_torque(enzymes_df.iloc[l]['superhelical']))
            if name == 'EXT_L' or name == 'EXT_R':
                s.append(0)
                c.append('white')
                m.append(MarkerStyle('o').get_path().transformed(MarkerStyle('o').get_transform()))
            #                m.append('o')
            else:
                s.append(enzyme_sizes[name])
                c.append(enzyme_colors[name])
                m.append(MarkerStyle(enzyme_shapes[name]).get_path().transformed(
                    MarkerStyle(enzyme_shapes[name]).get_transform()))
        #                m.append(enzyme_shapes[name])
        xl.append(x)
        yl.append(y)
        sl.append(s)
        cl.append(c)
        ml.append(m)
        sigma.append(sig)
        torque.append(tor)
        ttt = datetime.fromtimestamp(k * my_circuit.dt - 3600).strftime('%H:%M:%S')
        my_time.append('time=' + ttt)

    scat = ax[0].scatter(xl[0], yl[0], s=sl[0], c=cl[0], marker="o", ec='black', zorder=4)  # This plots RNAPs and NAPs

    lines = [ax[1].plot([], [], c=sigma_colour, lw=sigma_lw)[0] for _ in range(100)]  # This plots supercoiling

    if plot_torques:
        lines2 = [ax[2].plot([], [], c=torque_colour, lw=sigma_lw)[0] for _ in range(100)]  # This plots torque

    time_text = ax[0].text(0.0, 1.1, '', transform=ax[0].transAxes, fontsize=slabel)

    # ------------------------------------------------------------
    # ANIMATION
    # ------------------------------------------------------------
    def animate(i):
        x = xl[i]
        y = yl[i]
        s = sl[i]
        c = cl[i]
        m = ml[i]
        sig = sigma[i]
        tor = torque[i]
        xy = np.zeros((len(x), 2))
        xy[:, 0] = x
        xy[:, 1] = y
        scat.set_color(c)
        scat.set_sizes(s)
        scat.set_offsets(xy)
        scat.set_edgecolor('black')
        scat.set_paths(m)
        #        scat.set_paths([plt.Path.unit_regular_polygon(m)])  # Change marker shape to a square ('4')

        # scat.set_marker('o')

        # supercoiling stuff
        n = len(x)
        for j in range(10):
            lines[j].set_data([1, 1.02], [-.75, -.75])
            lines[j].set_linewidth(.2)

        for j in range(n):
            if j < n - 1:
                lines[j].set_data([x[j], x[j + 1]], [sig[j], sig[j]])  # , color=heat, lw=15 )
            if j == n - 1:
                lines[j].set_data([x[j], my_circuit.size], [sig[j], sig[j]])  # , color=heat, lw=15 )

            lines[j].set_linewidth(sigma_lw)

        # torque stuff
        if plot_torques:
            for j in range(10):
                lines2[j].set_data([1, 1.02], [-.75, -.75])
                lines2[j].set_linewidth(.2)

            for j in range(n):
                if j < n - 1:
                    lines2[j].set_data([x[j], x[j + 1]], [tor[j], tor[j]])  # , color=heat, lw=15 )
                if j == n - 1:
                    lines2[j].set_data([x[j], my_circuit.size], [tor[j], tor[j]])  # , color=heat, lw=15 )

                lines2[j].set_linewidth(sigma_lw)

        time_text.set_text(my_time[i])

        if plot_torques:
            return lines, lines2, scat, time_text
        else:
            return lines, scat, time_text

    # ANIMATE
    # -----------------------------------
    ani = animation.FuncAnimation(
        fig, animate, interval=10, frames=my_circuit.frames)  # , blit=True, frames=200, repeat=True)

    # SAVE OR SHOW
    # -----------------------------------

    writervideo = animation.FFMpegWriter(fps=30)
    # writervideo = animation.FFMpegWriter()
    #writervideo = animation.PillowWriter(fps=60)
    ani.save(output_file, writer=writervideo)
    # ani.save(output, writer=writervideo)
    # ani.save(output_file, writer='ImageMagick', fps=30)
    # ani.save(output_file, writer='FFMpeg', fps=30)
    # ani.save(output_file, writer='HTMLwriter', fps=30)
    # ani.save(output_file, fps=30)
    # plt.show()


# This function uses matplotlib ArtistAnimation instead of FuncAnimation. It also uses predifined png files to
# represent enzymes

def create_animation_linear_artist(my_circuit, sites_df, enzymes_df, output, out_format,
                                   site_type=None, site_colours=None):
    output_file = output + out_format
    h = 1.5
    dh = 0.5
    gx = np.array([0, my_circuit.size])
    gy = np.array([h, h])

    fig, ax = plt.subplots(2, figsize=(width * 1.2, height * 1.2), gridspec_kw={'height_ratios': [1, 2], 'hspace': 0.2})

    # Sizes
    # -----------------------------------
    ax[0].set_xlim(- 100, my_circuit.size + 100)
    ax[0].set_ylim(0, 2)
    ax[1].set_xlim(- 100, my_circuit.size + 100)
    ax[1].set_ylim(sigma_a, sigma_b)

    # labels and all that
    # -----------------------------------
    ax[0].grid(axis='x', zorder=1)
    ax[0].tick_params(labelleft=False, bottom=False, top=False)
    ax[1].grid(True, zorder=1)
    ax[1].set_xlabel("DNA (bp)", fontsize=slabel)
    ax[1].set_ylabel(r"$\sigma$", fontsize=slabel)

    # -----------------------------------
    # draw DNA
    # -----------------------------------
    ax[0].plot(gx, gy, lw=DNA_lw, color=DNA_colour, zorder=2)

    # -----------------------------------
    # Now draw genes
    # -----------------------------------
    for i, site in enumerate(my_circuit.site_list):
        if site.site_type == 'EXT':
            continue
        if site_type is not None and site.site_type != site_type:  # Only site types
            continue
        if 'DNA_' in site.site_type:  # I'm trying to avoid bare DNA binding sites
            continue

        x1 = site.end
        x0 = site.start
        dx = x1 - x0
        name = site.name
        if site_colours is not None:
            arrow = mpatches.FancyArrowPatch((x0, h), (x1, h), mutation_scale=25, color=site_colours[name], zorder=3,
                                             lw=gene_lw)
        else:
            arrow = mpatches.FancyArrowPatch((x0, h), (x1, h), mutation_scale=25, color=gene_colour, zorder=3,
                                             lw=gene_lw)
        ax[0].add_patch(arrow)
        if x0 < x1:
            a = x0 + abs(dx / 2)
        else:
            a = x1 + abs(dx / 2)
        ax[0].text(a, h - dh, name, fontsize=object_text)

    # -----------------------------------
    # Draw containers of mRNA
    # -----------------------------------
    for i, site in enumerate(my_circuit.site_list):
        if site.site_type == 'EXT':
            continue
        if site_type is not None and site.site_type != site_type:  # Only site types
            continue
        if 'DNA_' in site.site_type:  # I'm trying to avoid bare DNA binding sites
            continue

        rectangle, lit1, lit2 = create_mRNA_container(my_circuit=my_circuit, site=site)
        ax[0].add_patch(rectangle)
        ax[0].add_patch(lit1)
        ax[0].add_patch(lit2)

    # -----------------------------------
    # Draw enzyme labels
    # -----------------------------------
    y0 = 0.95#1.5#-.5
    text_y_offset= 0.03
    x =  .3 #10.0
    dx = .1 #1000 #0.5
    soom = 0.1
    #aax = ax.flat[0]
    for environmental in my_circuit.environmental_list:
        print(environmental.name)

        # Load image
        image = enzyme_shapes_dir + '/' + environmental.name + '.png'
        image_path = image

        if os.path.exists(image_path):
            img = plt.imread(image_path)
            im = OffsetImage(img, zoom=soom)

            # Annotate
            enzyme_annotation = AnnotationBbox(im, (x, y0), frameon=False,
                                               xycoords='figure fraction')
            fig.add_artist(enzyme_annotation)

            # Add text below the enzyme image
            text = f"{environmental.name}"
            ax[0].text(x, y0 - text_y_offset, text, ha='center', va='top', transform=fig.transFigure, fontsize=slabel*.75)

        else:
            raise ValueError(f"PNG file does not exist: {image_path}")

        x+=dx
        soom+=0.01


    # Data preparation before animation
    # -----------------------------------

    # Calculate max number of transcripts:
    # -----------------------------------
    gene_names = []
    for i, site in enumerate(my_circuit.site_list):
        if site.site_type == 'gene':
            gene_names.append(site.name)
    max_mRNA = -1
    for name in gene_names:
        mask = sites_df['name'] == name
        gene_df = sites_df[mask]
        n_mRNA = gene_df['unbinding'].sum()
        if n_mRNA > max_mRNA:
            max_mRNA = n_mRNA

    # Prepare df
    # -----------------------------------
    mask = sites_df['type'] == 'circuit'
    n_enzymes_df = sites_df[mask]

    # -----------------------------------
    # THE ANIMATION
    # -----------------------------------

    l = -1
    animation_frames = []  # Here, we will append the plots on each frame
    for k in range(my_circuit.frames):
        n_enz = n_enzymes_df.iloc[k]['#enzymes'] + 2  # 2 for the EXT

        # Add enzymes, collect sigma and positions
        # -----------------------------------------------------------------------------------------
        positions = []
        sigma = []
        molecule_drawings = []
        for i in range(n_enz):
            l = l + 1

            positions.append(enzymes_df.iloc[l]['position'])
            sigma.append(enzymes_df.iloc[l]['superhelical'])

            name = enzymes_df.iloc[l]['name']
            # Skip if EXT
            if name == 'EXT_L' or name == 'EXT_R':
                continue

            # Add molecule drawings
            # -----------------------------------------------------------------------------------------
            x = enzymes_df.iloc[l]['position']
            y = h

            # Call create_enzyme_annotation
            enzyme_annotation = create_enzyme_annotation(x, y, name)
            molecule_drawings.append(ax[0].add_artist(enzyme_annotation))  # And append it

        # Add superhelical density drawing
        # -----------------------------------------------------------------------------------------
        superhelical_drawings = []
        for j in range(n_enz):

            x1 = positions[j]
            y1 = sigma[j]
            y2 = sigma[j]
            if j < n_enz - 1:
                x2 = positions[j + 1]
            if j == n_enz - 1:
                x2 = my_circuit.size
            superhelical_drawings.append(ax[1].plot([x1, x2], [y1, y2], c=sigma_colour, lw=sigma_lw)[0])

        # Containers drawings
        # -----------------------------------------------------------------------------------------
        container_drawings = []
        for i, site in enumerate(my_circuit.site_list):
            if site.site_type == 'EXT':
                continue
            if site_type is not None and site.site_type != site_type:  # Only site types
                continue
            if 'DNA_' in site.site_type:  # I'm trying to avoid bare DNA binding sites
                continue

            # Calculate number of mRNA
            mask = (sites_df['name'] == site.name) & (sites_df['frame'] <= k)
            gene_df = sites_df[mask]
            n_mRNA = gene_df['unbinding'].sum()

            if n_mRNA < 1: continue # Skip if there are no mRNA yet

            if site_colours is not None:
                rectangle, lit1, lit2 = create_mRNA_container(my_circuit=my_circuit, site=site,
                                                       site_color=site_colours[site.name],
                                                       n_mRNA=n_mRNA, max_mRNA=max_mRNA)
            else:
                rectangle, lit1, lit2 = create_mRNA_container(my_circuit=my_circuit, site=site, site_color=gene_colour,
                                                       n_mRNA=n_mRNA, max_mRNA=max_mRNA)

            container_drawings.append(ax[0].add_patch(rectangle))
            container_drawings.append(ax[0].add_patch(lit1))
            container_drawings.append(ax[0].add_patch(lit2))

        # Time drawing
        # -----------------------------------------------------------------------------------------
        time_text = datetime.fromtimestamp(k * my_circuit.dt - 3600).strftime('%H:%M:%S')
        #        time_drawing = [ax[0].text(0.0, 1.1, time_text, transform=ax[0].transAxes, fontsize=slabel)[0]]
        time_drawing = ax[0].text(0.0, 1.1, time_text, transform=ax[0].transAxes, fontsize=slabel)

        # Join all animations
        # -----------------------------------------------------------------------------------------
        animation_frames.append(molecule_drawings + superhelical_drawings + container_drawings + [time_drawing])

    # ------------------------------------------------------------
    # ANIMATE
    # ------------------------------------------------------------
    ani = animation.ArtistAnimation(fig, animation_frames, interval=200, blit=True)

    # SAVE
    # -----------------------------------
    #writervideo = animation.FFMpegWriter(fps=30)
    #ani.save(output_file, writer=writervideo)

    # Let's try to make it faster
    # Using 'ultrafast' preset for faster encoding
    writervideo = animation.FFMpegWriter(fps=30, codec='libx264', extra_args=['-preset', 'ultrafast'])
    ani.save(output_file, writer=writervideo)


# Given an enzyme_name, it attaches a png representation to an ax object given the x and y positions.
# Note that the enzyme_name png representation (file) needs to be stored in effector_path
def create_enzyme_annotation(x, y, enzyme_name):  # , transform=None):
    image = enzyme_shapes_dir + '/' + enzyme_name + '.png'
    image_path = image

    soom = 0.1
    # Calculate the 5% range of the original number
    variation = 0.05 * soom
    # Generate a random number between -5% and +5%
    random_variation = random.uniform(0, variation)
    # Add the random variation to the original number
    soom += random_variation

    #if os.path.exists(image):
    #    im = OffsetImage(image, zoom=0.1)
    #else:
    #    raise ValueError('Effect png file does not exist:', image)
    #enzyme_annotation = AnnotationBbox(im, (x, y), frameon=False)

    if os.path.exists(image_path):
        img = plt.imread(image_path)
        im = OffsetImage(img, zoom=soom)
        #if transform is not None:
        #    enzyme_annotation = AnnotationBbox(im, (x, y), frameon=False, transform=transform)
        #else:
        enzyme_annotation = AnnotationBbox(im, (x, y), frameon=False)
        # ax.add_artist(enzyme_annotation)
    else:
        raise ValueError(f"PNG file does not exist: {image_path}")
    return enzyme_annotation


# Adds cartoon of container, which fills according the number of n_mRNA relative to a maximum number of max_mRNA
# Returns the rectangle and lit, which you'll need to add it outside the function
def create_mRNA_container(my_circuit, site, site_color=None, n_mRNA=None, max_mRNA=None):
    glass_w = my_circuit.size / 20.
    glass_h = .75
    glass_dup = .175
    glass_dho = glass_w
    y0 = 0.1

    # Get start position of glass
    if site.start < site.end:
        x0 = site.start + (site.end - site.start) / 2.
    else:
        x0 = site.end + (site.start - site.end) / 2.

    # Add empty glass; 3 lines + elipse
    # Lines
    #ax.plot([x0, x0], [glass_h, y0], '-k')
    #ax.plot([x0, x0 + glass_w], [y0, y0], '-k')
    #ax.plot([x0 + glass_w, x0 + glass_w], [y0, glass_h], '-k')
    # Lit
    lit_x = x0 + glass_dho / 2.

    if max_mRNA is not None and n_mRNA is not None:
        lit_h = glass_h * (n_mRNA / max_mRNA)
        rec_h = glass_h * (n_mRNA / max_mRNA) - glass_dup / 2.
    else:  # No mRNA yet
        lit_h = glass_h
        rec_h = glass_h - glass_dup / 2.

    if site_color is None:
        rectangle = mpatches.Rectangle([x0, y0], width=glass_w, height=rec_h, color=None, facecolor='white',
                                       edgecolor='black')
        lit1 = mpatches.Ellipse([lit_x, lit_h], width=glass_dho, height=glass_dup, color=None, facecolor='white',
                               edgecolor='black')
        lit2 = mpatches.Ellipse([lit_x, y0], width=glass_dho, height=glass_dup, color=None, facecolor='white',
                                edgecolor='black')
    else:
        rectangle = mpatches.Rectangle([x0, y0], width=glass_w, height=rec_h, facecolor=site_color,
                                       edgecolor='black', alpha=0.5)
        lit1 = mpatches.Ellipse([lit_x, lit_h], width=glass_dho, height=glass_dup,
                               facecolor=site_color, edgecolor='black', alpha=0.5)
        lit2 = mpatches.Ellipse([lit_x, y0], width=glass_dho, height=glass_dup,
                            facecolor=site_color, edgecolor='black', alpha=0.5)
    return rectangle, lit1, lit2
    #ax.add_patch(rectangle)
    #ax.add_patch(lit)
