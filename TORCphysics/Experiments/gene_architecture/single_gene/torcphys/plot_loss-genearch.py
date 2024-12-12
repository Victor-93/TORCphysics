import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
model_code='sus_GB-Stages-avgx2-02-'
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0
percentage_threshold = .1


loss_files = []

for pcase in promoter_cases:
    loss_files.append('susceptibility/' + model_code + pcase + '_dt' + str(dt) + '-values.csv')

#loss_files.append('calibrate_inferred-rates/GB-Stages-weak-kw0.02_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/GB-Stages-medium-kw0.02_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/GB-Stages-strong-kw0.02_dt1.0-values.csv')

#loss_files.append('calibrate_inferred-rates/Stages-weak_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/Stages-medium_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/Stages-strong_dt1.0-values.csv')

#loss_files.append('calibrate_inferred-rates/avgx2_simple_02/GB-Stages-avgx2-weak-kw0.02_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/avgx2_simple_02/GB-Stages-avgx2-medium-kw0.02_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/avgx2_simple_02/GB-Stages-avgx2-strong-kw0.02_dt1.0-values.csv')

# Erase this -----------------
#loss_files[0] = 'calibrate_inferred-rates/avgx2_simple_01/'+'GB-Stages-avgx2-' +'weak' + '-kw0.005_dt' + str(dt)+'-values.csv'
#loss_files[1] = 'calibrate_inferred-rates/avgx2_simple_01/'+'GB-Stages-avgx2-' + 'weak' + '-kw0.01_dt' + str(dt)+'-values.csv'
#loss_files[2] = 'calibrate_inferred-rates/avgx2_simple_01/'+'GB-Stages-avgx2-' + 'weak' + '-kw0.015_dt' + str(dt)+'-values.csv'

titles = ['Weak dt=1.0', 'Medium dt=1.0', 'Strong dt=1.0']

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
ms=5

# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we load
fig, axs = plt.subplots(len(loss_files), figsize=(width, len(loss_files)*height), tight_layout=True, sharex=True)

for i, loss_file in enumerate(loss_files):

    promoter = promoter_cases[i]

    if len(loss_files) >1:
        ax = axs[i]
    else:
        ax =axs

    df = pd.read_csv(loss_file)

    ax.set_title(titles[i])
#    ax.plot(df['test'], df['loss'], '-o', ms=ms, color='blue')

    df = df.sort_values(by='loss', ascending=False)  # , inplace=True)
    n = len(df['loss'])
    nconsidered = int(n * percentage_threshold)
    err_threshold = df['loss'].iloc[-nconsidered]
    print('Number of tests', n)
    print('Considered', nconsidered)
    print('For ', percentage_threshold * 100, '%')

    # Filter according error
    filtered_df = df[df['loss'] <= err_threshold]

    loss = df['loss'].to_numpy()
    floss = filtered_df['loss'].to_numpy()
    ax.plot(df['test'], df['loss'], 'o', ms=ms, color='blue', label='all')
    ax.plot(filtered_df['test'], filtered_df['loss'], 'o', ms=ms, color='red', label='best')

    ## Let's get the 5 smallest values:
    #df.sort_values(by=['loss'], ascending=False,inplace=True)
    #mdf = df.nsmallest(5, 'loss')
    #ax.plot(mdf['test'], mdf['loss'], 'o', ms=ms*1.5, color='red')
    #ax.plot(df['test'].iloc[0], mdf['loss'].iloc[0], 'o', ms=ms*2, color='green') # And the best

    #lo = df['loss'].to_numpy()
    #ax.plot( lo, 'o', ms=ms, color='blue')

    ax.grid(True)
    ax.set_xlabel('test')
    ax.set_ylabel('loss')
    ax.set_yscale('log')

    # Let's print some info:
    print(titles[i])
#    print(mdf)

    # Let's calculate averages and save the info:
    #-----------------------------------------------------------------------------------------------------------------------
    # We will store the dataframes in these lists so we can combine them at the end to create a table.
 #   avg_list = []
#    std_list = []

    # Drop the loss and test columns
    filtered_df = filtered_df.drop(columns=['loss', 'test'])

    # Calculate averages and standard deviations
    df_avg = filtered_df.mean(axis=0).to_frame().T.rename(index={0:'avg'})
    df_std = filtered_df.std(axis=0).to_frame().T.rename(index={0:'std'})

    # Save averages so we can load them for running
    df_avg.to_csv('gene-avg_'+model_code+ promoter+ '_dt' + str(dt) + '.csv', index=False, sep=',')

    # Add them to a list
#    avg_list.append(df_avg)
 #   std_list.append(df_std)

    # Join them for the table
#    df_avg = pd.concat([avg_list[0], avg_list[1]], axis=1)
#    df_std = pd.concat([std_list[0], std_list[1]], axis=1)
    new_df = pd.concat([df_avg, df_std], axis=0)
    new_df.to_csv('table_'+model_code+promoter+'_dt'+str(dt)+'.csv', index=False, sep=',')
#new_df.to_csv('table_dt'+str(dt)+'.csv', index=False, sep=',')


plt.show()
