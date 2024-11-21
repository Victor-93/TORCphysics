import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# This script does the pre-process that Aalap does, but only prints means for the yellow flourescent experiments.

# Function to perform pairwise t-tests
def perform_t_tests(df, variable):
    results = []
    groups = df.groupby('promoter')

    for promoter, group in groups:
        strains = group['strain'].unique()
        for i, strain1 in enumerate(strains):
            for strain2 in strains[i + 1:]:
                data1 = group[group['strain'] == strain1][variable]
                data2 = group[group['strain'] == strain2][variable]
                stat, pval = ttest_ind(data1, data2, nan_policy='omit')
                results.append({'promoter': promoter, 'strain1': strain1, 'strain2': strain2, 'p_value': pval})

    return pd.DataFrame(results)

# Import data
data = pd.read_csv("M_MT_MIT_S_ST_SIT_pTORC1_49_50_63_64_65_overnight_culture_24h_measurements.tsv", sep='\t')

# Group by 'date' and subtract the mean of the blanks
data['A600'] = data.groupby('date')['A600'].transform(lambda x: x - x[data['code'] == 'blank'].mean())
data['A700'] = data.groupby('date')['A700'].transform(lambda x: x - x[data['code'] == 'blank'].mean())
data['mhYFP'] = data.groupby('date')['mhYFP'].transform(lambda x: x - x[data['code'] == 'blank'].mean())
data['mRaspberry'] = data.groupby('date')['mRaspberry'].transform(lambda x: x - x[data['code'] == 'blank'].mean())

data['mhYFP_by_A600'] = data['mhYFP'] / data['A600']
data['mRas_by_A600'] = data['mRaspberry'] / data['A600']

data = data[data['code'] != 'blank']
data = data[data['code'] != 'pTORC66']

replicates = data.groupby(['bacterium', 'strain', 'plasmid']).size().reset_index(name='counts')
print(replicates)

#Analysis

# Filter E. coli data
data_mg1655 = data[data['bacterium'] == 'Escherichia coli K12 MG1655']

# Modify 'strain' column
data_mg1655['strain'] = data_mg1655['strain'].str.replace('FRT ', 'FRT\n')

# Define factor orders (can use categorical types in pandas)
code_order = ['pTORC1', 'pTORC49', 'pTORC50', 'pTORC65', 'pTORC63', 'pTORC64']
data_mg1655['code'] = pd.Categorical(data_mg1655['code'], categories=code_order, ordered=True)

promoter_order = ['None', 'PleuWT.1min mhYFP', 'PleuWT.1minRBS mhYFP', 'PleuWT.1 mhYFP', 'PleuWT.1min mRaspberry', 'PleuWT.1minRBS mRaspberry']
data_mg1655['promoter'] = pd.Categorical(data_mg1655['promoter'], categories=promoter_order, ordered=True)

strain_order = ['WT', 'ΔtopA::cat', 'ΔlacIZYA::FRT\nΔtopA::cat']
data_mg1655['strain'] = pd.Categorical(data_mg1655['strain'], categories=strain_order, ordered=True)

# Perform t-tests
t_tests_mhYFP = perform_t_tests(data_mg1655, 'mhYFP_by_A600')
t_tests_mRas = perform_t_tests(data_mg1655, 'mRas_by_A600')

# Adjust p-values using the Benjamini-Hochberg method
t_tests_mhYFP['p_adj'] = multipletests(t_tests_mhYFP['p_value'], method='fdr_bh')[1]
t_tests_mRas['p_adj'] = multipletests(t_tests_mRas['p_value'], method='fdr_bh')[1]

# Group and save so we can get the distributions!
# ----------------------------------------------------------------------------------------------------------------------
Ecoli_group = data_mg1655.groupby(['bacterium', 'promoter', 'strain']).apply(lambda x: x).reset_index(drop=True)
column_list = ['bacterium', 'promoter', 'strain', 'mhYFP_by_A600']
Ecoli_group = Ecoli_group[column_list]  # Select only the desired columns
Ecoli_group.to_csv('EColi_grouped.csv', index=False)

# Previous -------------------------------------------------------------------------------------------------------------
# Group by bacterium, promoter, and strain to get the mean
#mean_values_Ecoli = data_mg1655.groupby(['bacterium', 'promoter', 'strain'])['mhYFP_by_A600'].mean().reset_index()

#SAVE IT
#mean_values_Ecoli.to_csv('EColi_mhYFP_means.csv', index=False)

# New with std  --------------------------------------------------------------------------------------------------------
# Group by bacterium, promoter, and strain to get both the mean and standard deviation
mean_values_Ecoli = data_mg1655.groupby(['bacterium', 'promoter', 'strain'])['mhYFP_by_A600'].agg(
    mhYFP_by_A600='mean',  # Calculate mean
    mhYFP_by_A600_std='std'     # Calculate standard deviation
).reset_index()

# Save it to a CSV file
mean_values_Ecoli.to_csv('EColi_mhYFP_means.csv', index=False)

# SALMONELLA SECTION
# Filter data for 'Salmonella enterica Typhimurium SL1344'
data_sl1344 = data[data['bacterium'] == "Salmonella enterica Typhimurium SL1344"].copy()

# Modify strain names (replace "FRT " with "FRT\n" and "lacI" with "\nlacI")
data_sl1344['strain'] = data_sl1344['strain'].str.replace("FRT ", "FRT\n", regex=False)
data_sl1344['strain'] = data_sl1344['strain'].str.replace("lacI", "\nlacI", regex=False)

# Reordering 'code' factor levels
code_order = ["pTORC1", "pTORC49", "pTORC50", "pTORC65", "pTORC63", "pTORC64", "pTORC66"]
data_sl1344['code'] = pd.Categorical(data_sl1344['code'], categories=code_order, ordered=True)

# Reordering 'promoter' factor levels
promoter_order = ["None", "PleuWT.1min mhYFP", "PleuWT.1minRBS mhYFP", "PleuWT.1 mhYFP",
                  "PleuWT.1min mRaspberry", "PleuWT.1minRBS mRaspberry", "PleuWT.1 mRaspberry"]
data_sl1344['promoter'] = pd.Categorical(data_sl1344['promoter'], categories=promoter_order, ordered=True)

# Reordering 'strain' factor levels
strain_order = ["WT", "ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat", "ΔtopA::cat"]
data_sl1344['strain'] = pd.Categorical(data_sl1344['strain'], categories=strain_order, ordered=True)

# Group and save so we can get the distributions!
# ----------------------------------------------------------------------------------------------------------------------
Sal_group = data_sl1344.groupby(['bacterium', 'promoter', 'strain']).apply(lambda x: x).reset_index(drop=True)
Sal_group = Sal_group[column_list]  # Select only the desired columns
Sal_group.to_csv('Sal_grouped.csv', index=False)

# Previous -------------------------------------------------------------------------------------------------------------
# Group by bacterium, promoter, and strain to get the mean
#mean_values_Sal = data_sl1344.groupby(['bacterium', 'promoter', 'strain'])['mhYFP_by_A600'].mean().reset_index()

#SAVE IT
#mean_values_Sal.to_csv('Sal_mhYFP_means.csv', index=False)

# New with std  --------------------------------------------------------------------------------------------------------
# Group by bacterium, promoter, and strain to get both the mean and standard deviation
mean_values_Sal = data_sl1344.groupby(['bacterium', 'promoter', 'strain'])['mhYFP_by_A600'].agg(
    mhYFP_by_A600='mean',  # Calculate mean
    mhYFP_by_A600_std='std'     # Calculate standard deviation
).reset_index()

# Save it to a CSV file
mean_values_Sal.to_csv('Sal_mhYFP_means.csv', index=False)
