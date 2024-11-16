import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Description
#-----------------------------------------------------------------------------------------------------------------------
# We now want to plot the topoI, RNAP and gyrase KDEs for 2 or 3 key distances

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
#promoter_cases = ['weak', 'weak', 'strong']
promoter_cases = ['weak', 'medium', 'strong']
enzyme_names = ['topoI','RNAP']#, 'gyrase']
#dist_index = [0,9,11]
#dist_index = [0,10]
#dist_index = [10]
dist_index = [8,0]
time = 80000

dt=1.0#0.5

k_weak=0.02

n_inner_workers = 9
n_subsets = 2
n_sims = n_inner_workers * n_subsets  # This is the number of simulations launched to calculate unbinding rates

model_code = 'GB-Stages-'

experimental_files = []
calibration_files = []
for pcase in promoter_cases:
    # experimental_files.append('../../junier_data/inferred-rate_' + pcase + '.csv')
    experimental_files.append('../junier_data/inferred-rate_kw'+str(k_weak)+'_' + pcase + '.csv')
    #calibration_files.append('Stages-'+pcase+'_dt1.pkl')
    #calibration_files.append('Stages-'+pcase+'_dt'+str(dt)+'.pkl')
    #calibration_files.append('k_ini-Stages-'+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')
    # calibration_files.append('GB-Stages-'+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')

    # calibration_files.append(model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')
    calibration_files.append('calibrate_inferred-rates/reproduce-'+model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-o'
exp_ls = '--o'
titles = ['weak', 'medium', 'strong']

lines = ['-', '-', '-']
alphas = [1, 1, 1.0]
#enzyme_colors = {'RNAP': 'black', 'topoI': 'red', 'gyrase':'blue'}
enzyme_colors = {'RNAP': ['black', 'green'], 'topoI': ['red', 'orange'], 'gyrase':['blue', 'purple']}
colors = ['green', 'blue', 'red']

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
# Calculate rates as a function of distance
def get_prod_rates(results_list):
    x = [d['distance'] for d in results_list]
    # y = np.zeros_like(x)
    # ys = np.zeros_like(x)
    y = []
    ys = []
    # List with just results (not distance and not dict anymore)
    results_list2 = [d['result']['prod_rate'] for d in results_list]
    for j, rates_array in enumerate(results_list2): #case_results is the rate
        #rates = [d[0] for d in case_results]

        # Convert to a NumPy array
        #rates_array = np.array(case_re)

        # Calculate mean and standard deviation
        mean = rates_array[0]
        std = rates_array[1]
        #mean = np.mean(rates_array)
        #std = np.std(rates_array)
        y.append(mean)
        ys.append(std)

    rates = np.array([x, y, ys])
    return rates


# Organizes RNAP KDEs in a matrix format, thinking to plot it as a heatmap
def get_RNAP_KDE_matrix(results_list):
    # Let's define the matrix
    n = len(results_list)
    m = len(results_list[0]['result']['KDE']['RNAP']['kde_x'])
#    matrix = np.zeros((m, n))
    matrix = np.zeros((n, m))
    distance = []


    # Let's fill the matrix with KDE results
    for j in range(n):
        result_dict = results_list[j]
        distance.append(result_dict['distance'])
        kde_y = result_dict['result']['KDE']['RNAP']['kde_y']
        # matrix[:, j] = kde_y[0]  # 0-mean, 1-std
        matrix[j, :] = kde_y[0]  # 0-mean, 1-std

    max= 1# np.max(matrix)
    matrix = matrix/max
    return distance, matrix

# Organizes RNAP KDEs in lines format, thinking to plot it as a 3D collection
def get_RNAP_KDE_lines(results_list):
    # Let's define the matrix
    n = len(results_list)
    m = len(results_list[0]['result']['KDE']['RNAP']['kde_x'])
    # Will collect the lines
    RNAP_collection = [] #np.zeros((n, m))
    distance = []

    # Let's fill the matrix with KDE results
    for j in range(n):
        result_dict = results_list[j]
        distance.append(result_dict['distance'])
        y = result_dict['result']['KDE']['RNAP']['kde_y'][0]
        ys = result_dict['result']['KDE']['RNAP']['kde_y'][1]
        kde_x = result_dict['result']['KDE']['RNAP']['kde_x']
        x = kde_x-np.min(kde_x)

        RNAP_collection.append(np.array([x, y, ys]))
    return distance, RNAP_collection

# Get's the KDEs of RNAP, and topoisomerases for three key distances, and
# sets the start of the transcription unit (promoter) at x=0.
def get_enzymes_KDEs(results_list, names_list, index_list):
    output_list = []
    n = len(results_list)

    for i in index_list:
        out_dict = {}
        out_dict['index'] = i
        out_dict['distance'] = results_list[i]['distance']

        # We just need to obtain the value to align the enzymes
        x_disp =  np.min(results_list[i]['result']['KDE']['RNAP']['kde_x'])

        # Now we collect measurements
        for name in names_list:
            y = results_list[i]['result']['KDE'][name]['kde_y'][0]
            ys = results_list[i]['result']['KDE'][name]['kde_y'][1]
            kde_x = results_list[i]['result']['KDE'][name]['kde_x']
            x= kde_x - x_disp
            out_dict[name] = np.array([x, y, ys])
        output_list.append(out_dict)
    return output_list


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
    rates.append(get_prod_rates(data[0]['data']))

# Get KDEs
enzymes_KDEs = []
for i, data in enumerate(pickle_data):
    x=i
    enzymes_KDEs.append(get_enzymes_KDEs(data[0]['data'], enzyme_names, dist_index))

# Collect KDEs and arange them as heatmap
RNAP_KDEs = []
for i, data in enumerate(pickle_data):
    x=i
    RNAP_KDEs.append(get_RNAP_KDE_matrix(data[0]['data']))

# Collect KDEs and arange them as heatmap
RNAP_KDE_lines = []
for i, data in enumerate(pickle_data):
    x=i
    RNAP_KDE_lines.append(get_RNAP_KDE_lines(data[0]['data']))

# Plot - 3D lines
#-----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True, sharex=False)#, subplot_kw={'projection': '3d'})

for i in range(len(enzymes_KDEs)):

    ax = axs[i]
    axs[i].set_title(titles[i])

    # Get info
    case_kde = enzymes_KDEs[i]
    for j, kde_dict in  enumerate(case_kde):
        distance = kde_dict['distance']
        for name in enzyme_names:
            x = kde_dict[name][0]
            y = kde_dict[name][1]/time
            ys = kde_dict[name][2]/time#/1#/np.sqrt(time)
            if name == 'topoI' or name == 'gyrase':
                if j == 0:
                    x = x[15:-15]
                    y = y[15:-15]
                    ys = ys[15:-15]
                else:
                    x = x[:-15]
                    y = y[:-15]
                    ys = ys[:-15]

            label = name + ' ud:' +str(distance)
            #ax.plot(x, y, color=enzyme_colors[name], ls=lines[j], lw=lw,alpha=alphas[j])
            ax.plot(x, y, color=enzyme_colors[name][j], lw=lw, label=label)
            #ax.fill_between(x,y,y+ys,color=enzyme_colors[name][j], alpha=0.25)

    ax.set_ylabel(r'Density')
    ax.set_xlabel('Position (bp)')
    ax.grid(True)
#    ax.set_yscale('log')
    ax.legend(loc='best')
    #ax.set_ylim(0,100)

#fig.colorbar(im, )# ticks=[1, 2, 3])


plt.show()
#sys.exit()




