import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from TORCphysics import utils
from TORCphysics import binding_model as bm

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0#0.5

k_weak=0.003
# k_weak=0.0334

#model_code = 'GB-Stages-avg-'
model_code='sus_GB-Stages-avgx2-02-'
outfile = model_code+'fig'

#susceptibility_based = False
susceptibility_based = True

experimental_files = []
susceptibility_files = []
calibration_files = []
params_files = []
for pcase in promoter_cases:
    experimental_files.append('../junier_data/inferred-rate_kw'+str(k_weak)+'_' + pcase + '.csv')
    susceptibility_files.append('../junier_data/'+ pcase + '.csv')

    #calibration_files.append('calibrate_inferred-rates/avgx2_simple_02/'+model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')
    #params_files.append('calibrate_inferred-rates/avgx2_simple_02/'+model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.csv')

    # Susceptibility
    calibration_files.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    #calibration_files.append('susceptibility/reproduce-'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    #calibration_files.append('susceptibility/reproduce-'+model_code+pcase+'_dt'+str(dt)+'RNAPtracking_off.pkl')
    #RNAPtracking_off
    #params_files.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.csv')
    params_files.append('table_'+model_code+pcase+'_dt'+str(dt)+'.csv')
    #params_files.append('susceptibility/gene-avg_'+model_code+pcase+'_dt'+str(dt)+'.csv')

#calibration_files[1]='susceptibility/reproduce-'+model_code+promoter_cases[1]+'_dt'+str(dt)+'.pkl'

# Erase this -----------------
#experimental_files[0] = '../junier_data/inferred-rate_kw0.002_' + promoter_cases[0] + '.csv'
#experimental_files[1] = '../junier_data/inferred-rate_kw0.003_' + promoter_cases[1] + '.csv'
#experimental_files[2] = '../junier_data/inferred-rate_kw0.002_' + promoter_cases[2] + '.csv'

#calibration_files[0] = 'calibrate_inferred-rates/avgx2_simple_02/reproduce-' + model_code + promoter_cases[0] + '-kw0.002_dt' + str(dt) + '.pkl'
#calibration_files[1] = 'calibrate_inferred-rates/avgx2_simple_02/reproduce-' + model_code + promoter_cases[1] + '-kw0.003_dt' + str(dt) + '.pkl'
#calibration_files[2] = 'calibrate_inferred-rates/avgx2_simple_02/reproduce-' + model_code + promoter_cases[2] + '-kw0.002_dt' + str(dt) + '.pkl'

#params_files[0] = 'calibrate_inferred-rates/avgx2_simple_02/'+ model_code + promoter_cases[0] + '-kw0.002_dt' + str(dt)+'.csv'
#params_files[1] = 'calibrate_inferred-rates/avgx2_simple_02/'+ model_code + promoter_cases[1] + '-kw0.003_dt' + str(dt)+'.csv'
#params_files[2] = 'calibrate_inferred-rates/avgx2_simple_02/'+ model_code + promoter_cases[2] + '-kw0.002_dt' + str(dt)+'.csv'

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 7
height = 3.5
lw = 4
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-o'
exp_ls = '--o'
titles = ['Weak Promoter', 'Medium Promoter', 'Strong Promoter']

colors = ['green', 'blue', 'red']

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
# Calculate rates as a function of distance
def get_prod_rates(results_list):
    x = results_list['distance'] #[d['distance'] for d in results_list]
    y = []
    ys = []
    # List with just results (not distance and not dict anymore)
    for j, rates_array in enumerate(results_list['prod_rate']): #case_results is the rate

        # Calculate mean and standard deviation
        mean = rates_array[0]
        std = rates_array[1]
        y.append(mean)
        ys.append(std)

    rates = np.array([x, y, ys])
    return rates

# Calculate susceptibilities as a function of distance
def get_susceptibility(results_list):
    x = results_list['distance'] #[d['distance'] for d in results_list]
    y = []
    ys = []
    # List with just results (not distance and not dict anymore)
    for j, sus_array in enumerate(results_list['susceptibility']): #case_results is the rate

        # Calculate mean and standard deviation
        mean = sus_array[0]
        std = sus_array[1]
        y.append(mean)
        ys.append(std)

    rates = np.array([x, y, ys])
    return rates


def open_rate(superhelical, k_open, threshold, width):
    U = utils.opening_function(superhelical, threshold, width)
    rate = k_open * np.exp(-U)
    return rate

# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data = []
for pickle_file in calibration_files:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data.append(data)

# Calculate rates
rates = []
for i, data in enumerate(pickle_data):
    x=i
    rates.append(get_prod_rates(data[0]))

# Calculate susceptibilitie
susceptibility = []
for i, data in enumerate(pickle_data):
    x=i
    susceptibility.append(get_susceptibility(data[0]))

responses = []
for file in params_files:
    data = pd.read_csv(file)
    data = data.to_dict(orient='records')[0]  # Get the first (and only) dictionary
    responses.append(data)

# Plot
#-----------------------------------------------------------------------------------------------------------------------
sigma = np.arange(-.13, .0, 0.001)
#sigma = np.arange(-.2, .0, 0.001)
# Let's plot as we do the process
fig, axs = plt.subplots(3,2, figsize=(2*width, 3*height), tight_layout=True, sharex=False)

outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

if susceptibility_based:
    plot_array = susceptibility
    exp_files = susceptibility_files
else:
    plot_array = rates
    exp_files = experimental_files

for i, rate_array in enumerate(plot_array):
    ax = axs[i,0]

    ax.set_title(titles[i], fontsize=title_size) #color=colors[i],

    # Prepare calibration array and plot
    x = rate_array[0]
    y = rate_array[1]
    ys = rate_array[2]#/np.sqrt(n_sims)

    if i == 0:
        ax.plot(x, y, model_ls, lw=lw, color=colors[i], label='TORCPhysics')
    else:
        ax.plot(x, y, model_ls, lw=lw, color=colors[i])

    ax.fill_between(x, y-ys, y+ys, lw=lw, color=colors[i], alpha=0.2)

    # Load experimental and plot
    exp = pd.read_csv(exp_files[i]) # read
    x = exp['distance']
    y = exp['Signal']
    ys = exp['Error']
    if i == 0:
        ax.errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5,color=colors[i], label='Experiment')
    else:
        ax.errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5,color=colors[i])

    # Error quantification - meansquared error
    # error = np.sum((exp['Signal'] - rate_array[1]) ** 2)#/len(rate_array[1])
    error = np.sum((exp['Signal'] - rate_array[1]) ** 2)/len(rate_array[1])
    formatted_sum_err = "{:.3g}".format(error)  # Formatting to three significant figures
    er_label = r'$\epsilon={}$'.format(formatted_sum_err)

    # Add the text box to the plot
    props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
    ax.text(0.7, 0.1, er_label, transform=ax.transAxes, fontsize=font_size,
            verticalalignment='top', bbox=props)

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[i], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    #ax.set_ylabel(r'Expression Rate $(s^{-1})$', fontsize=font_size)
    ax.set_ylabel(r'Susceptibility', fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')
#    ax.set_ylim([0.75,1.4])



axs[0,0].legend(loc='best')

# Responses part
#-----------------------------------------------------------------------------------------------------------------------
x = sigma
promoter_label = ['Weak', 'Medium', 'Strong']
for i, resp_dict in enumerate(responses):

    # axs[i].set_title(titles[i])
    file = calibration_files[i]
    print(i,file)

    # Gaussian binding
    gy = bm.GaussianBinding(**resp_dict).rate_modulation(superhelical=x)
    y =  (gy - np.min(gy)) / (np.max(gy) - np.min(gy)) # Let's normalize
    axs[0,1].plot(x, y, lw=lw, color=colors[i])#, label=promoter_label[i])

    # Openning rate
    oy = open_rate(x, resp_dict['k_open'], resp_dict['threshold'], resp_dict['width'])
    y =  (oy - np.min(oy)) / (np.max(oy) - np.min(oy)) # Let's normalize
    #y=oy/resp_dict['k_open']
    axs[1,1].plot(x, y, lw=lw, color=colors[i], label=promoter_label[i])

# The spacer model
sb = bm.SpacerBinding(**{'k_on':1.0, 'superhelical_op': -0.06, 'spacer': 17})
y = sb.rate_modulation(superhelical=x)
axs[0,1].plot(x, y, lw=lw, color='gray', label='Spacer Model')

# Meyer model for PelE? model
mb= bm.MeyerPromoterOpening (**{'k_on':1.0})
y = mb.rate_modulation(superhelical=x)
y = (y - np.min(y)) / (np.max(y) - np.min(y))  # Let's normalize
axs[1,1].plot(x, y, lw=lw, color='gray', label='pelE')

axs[0,1].text(-0.1, 1.1, outside_label[3], transform=axs[0,1].transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')
axs[1,1].text(-0.1, 1.1, outside_label[4], transform=axs[1,1].transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

axs[0,1].grid(True)
axs[1,1].grid(True)
axs[0,1].legend(loc='best',fontsize=font_size)
axs[1,1].legend(loc='best',fontsize=font_size)
axs[0,1].set_xlabel('Superhelical Density', fontsize=xlabel_size)
axs[1,1].set_xlabel('Superhelical Density', fontsize=xlabel_size)
axs[0,1].set_ylabel(r'Rate Modulation', fontsize=xlabel_size)
axs[1,1].set_ylabel(r'Rate Modulation', fontsize=xlabel_size)
axs[0,1].set_title('Binding Rate (Gaussian Modulation)', fontsize=title_size)
axs[1,1].set_title('Open Complex Formation Rate (SIST)', fontsize=title_size)

# Bars of values
ax = axs[2,1]

nokeys = ['superhelical_op', 'spread', 'width', 'threshold', 'gamma', 'kappa', 'velocity', 'stall_torque']
# Remove the keys from all dictionaries in the list
for case_dict in responses:
    for key in nokeys:
        if key in case_dict:
            del case_dict[key]

keys = list(responses[0].keys())

# Number of cases and keys
n_cases = len(responses)
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
width = 0.8 / n_cases  # Dynamically calculate the width of each bar


# Plot each case in the list
for i, case_dict in enumerate(responses):

    values = [case_dict[key] for key in keys]  # Get values for the current case
    ax.bar(x + i * width - width * (n_cases - 1) / 2, values, width, label=promoter_cases[i], color=colors[i])

# Add labels, title, and custom ticks
ax.set_xlabel('Rate Type', fontsize=xlabel_size)
ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
ax.set_title('Transition Rates', fontsize=title_size)
ax.set_xticks(x)
ax.set_xticklabels(keys)

axs[2,1].text(-0.1, 1.1, outside_label[5], transform=axs[2,1].transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')


#plt.savefig(model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'_v2.png')
plt.savefig(outfile+'.pdf')

plt.show()



