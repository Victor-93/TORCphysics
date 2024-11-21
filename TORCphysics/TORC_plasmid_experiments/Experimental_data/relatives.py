# You have to load the means, and define a reference. With that, try plotting the bars.
# You have to choose which would be the reference system (the most basal one), and use this to compare everything else.
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#  Description
#-----------------------------------------------------------------------------------------------------------------------
# I want to produce reference curves, so I can use the relationship for the modelling using TORCPhysics.
# Here, I will use the Salmonella bacterium as the reference system.
# You can change the reference system (file), strain and promoter (dict)

# Input data
#-----------------------------------------------------------------------------------------------------------------------
Ecoli_input = 'EColi_mhYFP_means.csv'
Ecoli_input_group = 'EColi_grouped.csv'
Sal_input = 'Sal_mhYFP_means.csv'
Sal_input_group = 'Sal_grouped.csv'

reference_dict = {'file':Sal_input, 'promoter': 'PleuWT.1min mhYFP', 'strain': 'WT'}

# This is a list of the promoters I want to filter - The ones I'll end up using
promoter_target =  ['PleuWT.1min mhYFP', 'PleuWT.1 mhYFP']

out_file = 'reference'
# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 5
height = 5
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-'
exp_ls = '--'
titles = ['weak', 'medium', 'strong']

colors = ['green', 'blue', 'red']

# Process - Using averages
#-----------------------------------------------------------------------------------------------------------------------
# Load Data
Ecoli_data = pd.read_csv(Ecoli_input)
Sal_data = pd.read_csv(Sal_input)

# Get reference value and data
df = pd.read_csv(reference_dict['file']) # Load the reference file
#reference_value = df.loc[(df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'mhYFP_by_A600'].values[0]
reference_value = df.loc[(df['promoter'] == reference_dict['promoter']) & (df['strain'] == reference_dict['strain']), 'mhYFP_by_A600'].values[0]

# Let' transform the data using the reference
Ecoli_data['mhYFP_by_A600'] = Ecoli_data['mhYFP_by_A600']/reference_value
Ecoli_data['mhYFP_by_A600_std'] = Ecoli_data['mhYFP_by_A600_std']/reference_value
Sal_data['mhYFP_by_A600'] = Sal_data['mhYFP_by_A600']/reference_value
Sal_data['mhYFP_by_A600_std'] = Sal_data['mhYFP_by_A600_std']/reference_value

# Filter the DataFrame using the 'promoter' column and `isin()`
Ecoli_data = Ecoli_data.loc[Ecoli_data['promoter'].isin(promoter_target)]
Sal_data = Sal_data.loc[Sal_data['promoter'].isin(promoter_target)]

# And join them to form a new reference df
ref_df = pd.concat([Ecoli_data, Sal_data], ignore_index=True)
ref_df = ref_df.rename(columns={'mhYFP_by_A600': 'relative'}) # Rename the column with measurement
ref_df = ref_df.rename(columns={'mhYFP_by_A600_std': 'std'}) # Rename the column with measurement
ref_df.to_csv(out_file+'.csv', index=False)

# Process - Using distributions instead - (It is exactly the same process)
#-----------------------------------------------------------------------------------------------------------------------
# Load Data
Ecoli_data = pd.read_csv(Ecoli_input_group)
Sal_data = pd.read_csv(Sal_input_group)

# Let' transform the data using the reference
Ecoli_data['mhYFP_by_A600'] = Ecoli_data['mhYFP_by_A600']/reference_value
Sal_data['mhYFP_by_A600'] = Sal_data['mhYFP_by_A600']/reference_value

# Filter the DataFrame using the 'promoter' column and `isin()`
Ecoli_data = Ecoli_data.loc[Ecoli_data['promoter'].isin(promoter_target)]
Sal_data = Sal_data.loc[Sal_data['promoter'].isin(promoter_target)]

# And join them to form a new reference df
dist_ref_df = pd.concat([Ecoli_data, Sal_data], ignore_index=True)
dist_ref_df = dist_ref_df.rename(columns={'mhYFP_by_A600': 'relative'}) # Rename the column with measurement
dist_ref_df.to_csv(out_file+'_distribution.csv', index=False)

# And plot - for average
#-----------------------------------------------------------------------------------------------------------------------
# Create a FacetGrid to organize the rows by 'bacterium' and group bars by 'promoter'
df = ref_df
g = sns.catplot(
    data=df, kind='bar',
    x='strain', y='relative', hue='promoter',
    col='bacterium', col_wrap=2,  # One column per bacterium
    height=height, aspect=1
    #    palette='coolwarm'  # Using the 'coolwarm' colormap
)

# Then, add error bars manually
for ax in g.axes.flat:
    for i, bar in enumerate(ax.patches):
        # Add error bars to each bar using the std values
        row = df.iloc[i]  # Get the corresponding row in df
        std_val = row['std']  # Get the std value for that bar
        ax.errorbar(
            bar.get_x() + bar.get_width() / 2,  # x position of the bar's center
            bar.get_height(),  # y value at the top of the bar
            yerr=std_val,  # Error bar from the 'std' column
            fmt='none',  # No marker for the error bar itself
            ecolor='black',  # Color of the error bar
            capsize=5  # Size of the caps on the error bars
        )
        print(row)

# Customize the plot
g.set_axis_labels('Strain', 'Relative Value')
g.set_titles('Bacterium: {col_name}')
g.despine(left=True)

# Add grid lines to all subplots
for ax in g.axes.flat:
    ax.grid(True)  # Enable grid lines for each subplot

g.savefig(out_file+".png", dpi=300, bbox_inches='tight')  # Save as PNG with high DPI

# Show the plot
plt.show()

# And plot - for distribution
#-----------------------------------------------------------------------------------------------------------------------
# Create a FacetGrid to organize the rows by 'bacterium' and group bars by 'promoter'
df = dist_ref_df
g = sns.catplot(
    data=df, kind='bar',
    x='strain', y='relative', hue='promoter',
    col='bacterium', col_wrap=2,  # One column per bacterium
    height=height, aspect=1
    #    palette='coolwarm'  # Using the 'coolwarm' colormap
)


# Customize the plot
g.set_axis_labels('Strain', 'Relative Value')
g.set_titles('Bacterium: {col_name}')
g.despine(left=True)

# Add grid lines to all subplots
for ax in g.axes.flat:
    ax.grid(True)  # Enable grid lines for each subplot

g.savefig(out_file+"_dist.png", dpi=300, bbox_inches='tight')  # Save as PNG with high DPI

# Show the plot
plt.show()





