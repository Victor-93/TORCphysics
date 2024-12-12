import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Description
#-----------------------------------------------------------------------------------------------------------------------
# The porpuse of this script is to prepare the experimental curves for calibration.
# Rather than directly calibrating using susceptibilities (since they cause a lot of variation), we will directly
# use inferred expression rates. The idea is to use the experimental relative expression rate of the weak promoter
# as reference, assign it a expression rate, and from the scale it to the other two (medium & strong) experimental
# curves.
# The process is as follows. We choose a expression rate for the weak promoter k_weak. From the relative expression
# rates from Junier data, we know the relationship between the relative expression rates. Taking as the
# 5th data-point as reference (as they do it), the medium promoter is 2.772 times stronger than the weak promoter,
# and the strong promoter is 4.550 times stronger than the weak. With this in mean, we re-scale the susceptibility
# curves to expression rates as follows: f_weak = s_weak* k_weak; f_medium = s_medium * 2.772 * k_weak;
# f_strong = s_strong * 4.550 * k_weak, where s_* denotes the susceptibility curves, and f_* the resulting
# expression rate curves.

# Inputs
#-----------------------------------------------------------------------------------------------------------------------

# Junier experimental data
s_weak = 'weak.csv'
s_medium = 'medium.csv'
s_strong = 'strong.csv'

# Factors
weak_F = 1.0
medium_F = 2.772
strong_F = 4.550
factors = [weak_F, medium_F, strong_F]

# Weak promoter expression rate
#k_weak = 0.0334#0.02#0.02 # 0.02 ~ 2 transcripts per 100s, or 1 transcript every 50seconds
k_weak = 0.005#0.0275 #0.0225# 0.03 #0.025

# Output parms
out_prefix = 'inferred-rate_kw'+str(k_weak)+'_'

# Plotting params
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

# Process & plot
#-----------------------------------------------------------------------------------------------------------------------
cinput = [s_weak, s_medium, s_strong]

titles = ['weak', 'medium', 'strong']

colors = ['green', 'blue', 'red']

# Let's plot as we do the process
fig, axs = plt.subplots(4, figsize=(width, 4*height), tight_layout=True, sharex=True)

for i, cin in enumerate(cinput):

    # Load and process
    susceptibility = pd.read_csv(cin) # read
    coutput = out_prefix + cin
    rate = susceptibility.copy()
    rate['Signal'] = rate['Signal'] * k_weak * factors[i]
    rate['Error'] = rate['Error'] * k_weak * factors[i]
    rate.to_csv(coutput, index=False)  # output our new curves

    # And plot them to visualize them
    x = rate['distance']
    y = rate['Signal']
    ys = rate['Error']
    axs[i].set_title(titles[i])
    axs[i].plot(x, y, '-o', lw=lw, color=colors[i])
    #axs[i].fill_between(x, y-ys, y+ys, 'black')
    axs[i].set_xlabel('Distance')
    axs[i].set_ylabel('expression rate')
    axs[i].grid(True)
    axs[i].set_xscale('log')

    # Last one
    axs[3].plot(x, y, '-o', lw=lw, color=colors[i])

# Last one
axs[3].set_title('all together')
axs[3].set_xlabel('Distance')
axs[3].set_ylabel('expression rate')
axs[3].grid(True)
axs[3].set_xscale('log')

plt.show()


