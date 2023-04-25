import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import matplotlib.patches as mpatches
from TORCphysics import analysis as an

# This module is produced to make one's life easier and be able to quickly process results.
# With this module you can create animations,

# Figure parmameters
width = 10
height = 6

# colors
enzyme_colors = {'RNAP': 'blue', 'IHF': 'yellow', 'FIS': 'red', 'lacI': 'black', 'ori': 'yellow'}
gene_colour = 'green'
DNA_colour = 'black'
sigma_colour = 'red'
topoI_colour = 'red'
gyrase_colour = 'cyan'

# Sizes
enzyme_sizes = {'RNAP': 300, 'IHF': 500, 'FIS': 500, 'lacI': 500, 'ori': 500}
RNAP_size = 300
FIS_size = 500
IHF_size = 500
LAC_size = 500
ORI_size = 500
DNA_lw = 12
gene_lw = 5
sigma_lw = 5

# Shapes
enzyme_shapes = {'RNAP': 'o', 'IHF': 'o', 'FIS': 'o', 'lacI': 's', 'ori': 's'}

# text size
slabel = 15
object_text = 10  # NAPs, genes, etc...


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
# TODO: check how to use kwargs and all that
# TODO: Also add the ignore feature to other relevant functions
# TODO: Maybe one wants to include specific sites names? instead of ignoring?
# TODO: You need to test these functions
# TODO: Change ta,tb for fa,fb. Because ta is time, and we give it in frames

def ax_params(axis, xl, yl, grid, legend):
    axis.grid(grid)
    axis.set_ylabel(yl)
    axis.set_xlabel(xl)
    if legend:
        axis.legend(loc='best')


# Plots the steady state curve: log( sum of initiation events / time ) vs time
def plot_steady_state_initiation_curve(my_circuit, sites_df, axs, ignore=None,
                                       ta=0, tb=-1, colors=None, site_type=None, labels=True):
    time = np.arange(0, my_circuit.dt * (my_circuit.frames + 1), my_circuit.dt)
    if site_type is None:
        curves, rates, names = an.initiation_rates(sites_df, time, ta, tb)
    else:
        curves, rates, names = an.initiation_rates_by_type(sites_df, site_type, time, ta, tb)

    for i, name in enumerate(names):
        curve = curves[i]
        rate = rates[i]
        if ignore is not None:
            if name in ignore:
                continue
        my_label = name + f' k={rate:.4f}s'
        if colors is not None:
            axs.plot(time, curve, color=colors[name], label=my_label, alpha=0.5)
        else:
            axs.plot(time, curve, label=my_label, alpha=0.5)
    if labels:
        ax_params(axis=axs, xl='time (seconds)', yl=r'$\log$($\sum$initiation)/time', grid=True, legend=True)


# TODO: Cross-correlation with no specific site? That might be a mess
# Plot cross-correlations between sites with respect another site with name ref_name
# t0 = Time in which we assume the system has reached the steady state
# ref_name is the name of the site that you want to calculate the cross-correlation with. It has to have the same type
# as site_type in case you specify.
def plot_cross_correlation_with_site(my_circuit, sites_df, ref_name, axs, ignore=None,
                                     ta=0, tb=-1, colors=None, site_type=None, labels=True):
    if site_type is None:
        signals, names = an.build_signals(sites_df)
    else:
        signals, names = an.build_signal_by_type(sites_df, site_type)
    time = np.arange(0, my_circuit.dt * len(signals[0]), my_circuit.dt)
    signals_t0 = []
    for signal in signals:
        signals_t0.append(signal[ta:tb])
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
            axs.plot(lag, cross[index, i, :], color=colors[name], label=my_label)
        else:
            axs.plot(lag, cross[index, i, :], label=my_label)
    if labels:
        ax_params(axis=axs, xl='time lag (seconds)', yl='Cross-correlation w ' + ref_name, grid=True, legend=True)
    axs.set_xlim(-200, 200)


# Plots supercoiling profiles at the sites and global
def plot_supercoiling_profiles(my_circuit, sites_df, axs, colors=None, site_type=None, labels=True):
    time = np.arange(0, my_circuit.dt * (my_circuit.frames + 1), my_circuit.dt)
    # Let's plot the global superhelical density
    mask = sites_df['type'] == 'circuit'  # This one contains global superhelical density
    global_sigma = sites_df[mask]['superhelical'].to_numpy()
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
        mask = sites_df['name'] == name
        superhelical = sites_df[mask]['superhelical'].to_numpy()
        if colors is not None:
            axs.plot(time, superhelical, color=colors[name], label=name, alpha=0.5)
        else:
            axs.plot(time, superhelical, label=name, alpha=0.5)
    axs.plot(time, global_sigma, color='black', label='global')  # and the global

    ax_params(axis=axs, xl='time (seconds)', yl='Supercoiling at site', grid=True, legend=True)


# This one plots the signal profiles.
def plot_signal_profiles(my_circuit, sites_df, axs, colors=None, site_type=None, labels=True):
    if site_type is None:
        signals, names = an.build_signals(sites_df)
    else:
        signals, names = an.build_signal_by_type(sites_df, site_type)
    time = np.arange(0, my_circuit.dt * len(signals[0]), my_circuit.dt)
    for i, signal in enumerate(signals):
        name = names[i]
        if colors is not None:
            axs.plot(time, signal, color=colors[name], label=names[i], alpha=0.5)
        else:
            axs.plot(time, signal, label=names[i], alpha=0.5)
    if labels:
        ax_params(axis=axs, xl='time (seconds)', yl='Transcription signal', grid=True, legend=True)


# Sites rate curves - Let's plot the rates modulated by supercoiling
def plot_site_response_curves(my_circuit, axs, colors=None, site_type=None, labels=True):
    i = -1
    for site in my_circuit.site_list:
        if site_type is None or site.site_type == site_type:
            i += 1
            rate, x = an.site_activity_curves(site)
            if colors is not None:
                axs.plot(x, rate, color=colors[site.name], label=site.name)
            else:
                axs.plot(x, rate, label=site.name)

    if labels:
        ax_params(axis=axs, xl=r'$\sigma$', yl=r'Initiation rate ($s^{-1}$)', grid=True, legend=True)


# Plots the topoisomerase activity curves of a continuum model
def plot_topoisomerase_activity_curves_continuum(my_circuit, axs, labels=True):
    i = -1
    for environmental in my_circuit.environmental_list:
        if environmental.enzyme_type == 'topo':
            i += 1
            topo_curve, x = an.topoisomerase_activity_curves_continuum(environmental, dt=my_circuit.dt)
            axs.plot(x, topo_curve, label=environmental.name)
            if i == 0:
                topo_sum = np.zeros_like(topo_curve)
            topo_sum += topo_curve
    axs.plot(x, topo_sum, color='black', label='sum')
    if labels:
        ax_params(axis=axs, xl=r'$\sigma$', yl=r'$\sigma$ removed per timestep', grid=True, legend=True)


def create_animation_linear(my_circuit, sites_df, enzymes_df, output, out_format,
                            site_type=None, site_colours=None):
    output_file = output + out_format
    h = 1.5
    dh = 0.5
    gx = np.array([0, my_circuit.size])
    gy = np.array([h, h])

    fig, ax = plt.subplots(2, figsize=(width, height), gridspec_kw={'height_ratios': [1, 2], 'hspace': 0.2})

    # Sizes
    # -----------------------------------
    ax[0].set_xlim(- 100, my_circuit.size + 100)
    ax[0].set_ylim(0, 2)
    ax[1].set_xlim(- 100, my_circuit.size + 100)
    ax[1].set_ylim(-0.25, .25)

    # labels and all that
    # -----------------------------------
    ax[0].grid(axis='x', zorder=1)
    ax[1].grid(True, zorder=1)
    # ax[0].set_xlabel("DNA (bp)", fontsize=slabel)
    ax[1].set_xlabel("DNA (bp)", fontsize=slabel)
    ax[1].set_ylabel(r"$\sigma$", fontsize=slabel)
    ax[0].tick_params(labelleft=False, bottom=False, top=False)

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
        for i in range(n_enz):
            l = l + 1
            name = enzymes_df.iloc[l]['name']
            x.append(enzymes_df.iloc[l]['position'])
            y.append(h)
            sig.append(enzymes_df.iloc[l]['superhelical'])
            if name == 'EXT_L' or name == 'EXT_R':
                s.append(0)
                c.append('white')
                m.append('o')
            else:
                s.append(enzyme_sizes[name])
                c.append(enzyme_colors[name])
                m.append(enzyme_shapes[name])
        xl.append(x)
        yl.append(y)
        sl.append(s)
        cl.append(c)
        ml.append(m)
        sigma.append(sig)

    scat = ax[0].scatter(xl[0], yl[0], s=sl[0], c=cl[0], marker="o", zorder=4)  # This plots RNAPs and NAPs

    lines = [ax[1].plot([], [], c=sigma_colour, lw=sigma_lw)[0] for _ in range(10)]  # This plots supercoiling

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
        xy = np.zeros((len(x), 2))
        xy[:, 0] = x
        xy[:, 1] = y
        scat.set_color(c)
        scat.set_sizes(s)
        scat.set_offsets(xy)

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

        return lines, scat

    # ANIMATE
    # -----------------------------------
    ani = animation.FuncAnimation(
        fig, animate, interval=10, frames=my_circuit.frames)  # , blit=True, frames=200, repeat=True)

    # SAVE OR SHOW
    # -----------------------------------
    #writervideo = animation.FFMpegWriter(fps=60)
    #ani.save(output_file, writer=writervideo)
    #ani.save(output_file, writer='imagemagick', fps=30)
    ani.save(output_file, writer='Pillow', fps=30)
    #ani.save(output_file, writer='HTMLwriter', fps=30)
    #ani.save(output_file, fps=30)
    #plt.show()

