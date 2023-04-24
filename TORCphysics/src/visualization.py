import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import matplotlib.patches as mpatches

# This module is produced to make one's life easier and be able to quickly process results.
# With this module you can create animations,

# Figure parmameters
width = 10
height = 6

# colours
RNAP_colour = 'blue'
IHF_colour = 'purple'
FIS_colour = 'red'
gene_colour = 'green'
DNA_colour = 'black'
LAC_colour = 'gray'
ORI_colour = 'yellow'
sigma_colour = 'red'

# Sizes
RNAP_size = 300
FIS_size = 500
IHF_size = 500
LAC_size = 500
ORI_size = 500
DNA_lw = 12
gene_lw = 5
sigma_lw = 5

# text size
slabel = 15
object_text = 10  # NAPs, genes, etc...

def create_animation_linear(my_circuit, sites_df, enzymes_df, frames, output, out_format):
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
    n_genes = len(my_circuit.site_list) - 2
    for i in range(n_genes):
        x1 = my_circuit.site_list[i + 2].end
        x0 = my_circuit.site_list[i + 2].start
        dx = x1 - x0
        name = my_circuit.site_list[i + 2].name
        arrow = mpatches.FancyArrowPatch((x0, h), (x1, h), mutation_scale=25, color=gene_colour, zorder=3, lw=gene_lw)
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
    for k in range(frames):
        n_enz = n_enzymes_df.iloc[k]['#enzymes'] + 2  # 2 for the EXT
        x = []
        y = []
        s = []
        c = []
        m = []
        sig = []
        for i in range(n_enz):
            l = l + 1
            x.append(enzymes_df.iloc[l]['position'])
            y.append(h)
            sig.append(enzymes_df.iloc[l]['superhelical'])
            name = enzymes_df.iloc[l]['name']
            if name == 'RNAP':
                s.append(RNAP_size)
                c.append(RNAP_colour)
                m.append("o")
            if name == 'IHF':
                s.append(IHF_size)
                c.append(IHF_colour)
                m.append("s")
            if name == 'FIS':
                s.append(FIS_size)
                c.append(FIS_colour)
                m.append("o")
            if name == 'EXT_L' or name == 'EXT_R': #EXTRA
                s.append(0)
                c.append('black')
                m.append("o")
            if name == 'lacI':
                s.append(LAC_size)
                c.append(LAC_colour)
                m.append("s")
            if name == 'ori':
                s.append(ORI_size)
                c.append(ORI_colour)
                m.append("s")
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
        fig, animate, interval=10, frames=frames)  # , blit=True, frames=200, repeat=True)

    # SAVE OR SHOW
    # -----------------------------------
#    ani.save(output_file, writer='imagemagick', fps=30)
    ani.save(output_file, writer='Pillow', fps=30)

