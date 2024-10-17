import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
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

# TODO: Let's first plot the 3D densities for RNAPs. Then we might do it for topoI and gyrase
#  We alos need to compare the parameters
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
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True, sharex=True)#, subplot_kw={'projection': '3d'})

for i in range(len(RNAP_KDE_lines)):

    ax = axs[i]
    axs[i].set_title(titles[i])

    distance = RNAP_KDE_lines[i][0]
    RNAP_collection = RNAP_KDE_lines[i][1]

    # Plot each 3D plot in separate subplots
    for j, dist in enumerate(distance):
        x = RNAP_collection[j][0]
        y = RNAP_collection[j][1]
        ys = RNAP_collection[j][2]
        # ax.plot(x, y, zs=dist)
        ax.plot(x, y)

    #ax.set_ylabel(r'Upstream distance')
    #ax.set_xlabel('Transcription unit')
    ax.grid(True)

#fig.colorbar(im, )# ticks=[1, 2, 3])


plt.show()
#sys.exit()
# Plot matrix
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we do the process
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True, sharex=True)

gene_length = 900
x_pos =np.array([0,25,50,75,100])
x_labels = gene_length*x_pos/100

# Let's find the maximum
mymax=0
for i in range(len(RNAP_KDEs)):
    RNAP_matrix = RNAP_KDEs[i][1]

    max = np.max(RNAP_matrix)
    if max > mymax:
        mymax = max

for i in range(len(RNAP_KDEs)):

    ax = axs[i]
    axs[i].set_title(titles[i])

    distance = RNAP_KDEs[i][0]
    RNAP_matrix = RNAP_KDEs[i][1]/mymax

    im  = ax.imshow( RNAP_matrix, cmap= 'viridis', interpolation= 'nearest', aspect=5, vmin=0, vmax=1)

    y_pos = np.arange(len(distance))
    y_labels = distance

    ax.set_yticks( y_pos )
    ax.set_yticklabels( y_labels )

    ax.set_xticks( x_pos )
    ax.set_xticklabels( x_labels )

    ax.set_ylabel(r'Upstream distance')
    ax.set_xlabel('Transcription unit')
    ax.grid(True)

#fig.colorbar(im, )# ticks=[1, 2, 3])


plt.show()



