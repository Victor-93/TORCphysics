import sys
#sys.path.append("/users/ph1vvb")
from TORCphysics import params
from TORCphysics import parameter_search as ps
import pandas as pd
import numpy as np
import pickle

# Description
# ----------------------------------------------------------------------------------------------------------------------
# I want to automatically create oparams and circuit files for TORCPhysics following the optimization process.

# Model & simulation conditions
# ----------------------------------------------------------------------------------------------------------------------
lacI_blocking = True

file_in = 'op_TORC_plasmid.csv'

if lacI_blocking == True:
    file_in = 'block-'+file_in

gene_names = ['PleuWT', 'tetA', 'bla', 'antitet']
environment_filename = 'environment.csv'

# Read inputs
# ----------------------------------------------------------------------------------------------------------------------
df_oparams = pd.read_csv(file_in)
df_environment = pd.read_csv(environment_filename)

# We need to create environments for Ecoli and Salmonella, and we include RNAP, topos and lacI

# Create oparams
# ----------------------------------------------------------------------------------------------------------------------
# For genes
for gene in gene_names:
    d = df_oparams.filter(like=gene) # responses from optimization
    s = pd.read_csv('../responses/'+gene+'.csv') # Responses from SIST

    my_oparam = pd.DataFrame({
        'k_on': d[gene + '_k_on'], 'k_off': d[gene + '_k_off'].iloc[0],
        'k_open': d[gene + '_k_open'].iloc[0], 'k_closed': d[gene + '_k_closed'].iloc[0],
        'k_ini': d[gene + '_k_ini'].iloc[0], 'superhelical_op': d[gene + '_superhelical_op'].iloc[0],
        'spread': d[gene + '_spread'].iloc[0],
        'width': s['width'].iloc[0], 'threshold': s['threshold'].iloc[0]
    })
    # And let's save it
    my_oparam.to_csv('oparam_'+gene+'.csv', index=False,sep=',')

# For lacI
my_oparam = pd.DataFrame({
    'k_on': df_oparams['lacI_k_on'],
    'k_off': df_oparams['lacI_k_off'],
    'k_bridge_on': df_oparams['bridge_on'],
    'k_bridge_off': df_oparams['bridge_off'],
    'leakage': df_oparams['leakage']
})
my_oparam.to_csv('oparam_lacI.csv', index=False, sep=',')

# TODO: Create environment
# ----------------------------------------------------------------------------------------------------------------------
# Let's use the df_environment as template

# Apply lacI oparam location
lacI_mask = df_environment['name'] == 'lacI'
df_environment.loc[lacI_mask, 'unbinding_model'] = 'LacISimpleUnBinding'
df_environment.loc[lacI_mask, 'unbinding_oparams'] = '../optimization/oparam_lacI.csv'
df_environment.loc[lacI_mask, 'effect_model'] = 'LacIPoissonBridging'
df_environment.loc[lacI_mask, 'effect_oparams'] = '../optimization/oparam_lacI.csv'

# Let's apply the model for RNAP
RNAP_mask = df_environment['name'] == 'RNAP'
df_environment.loc[RNAP_mask, 'unbinding_model'] = 'RNAPStagesSimpleUnbindingv2'
df_environment.loc[RNAP_mask, 'effect_model'] = 'RNAPStagesStallv2'

# E coli
# --------------------------------------------
df_environment.to_csv('Ecoli_environment.csv', index=False, sep=',')

# Salmonella
# --------------------------------------------
# Change concentration of topoisomerases
mask = df_environment['name'] == 'topoI'
df_environment.loc[mask, 'concentration'] = df_oparams['topoI_concentration'].iloc[0]
mask = df_environment['name'] == 'gyrase'
df_environment.loc[mask, 'concentration'] = df_oparams['gyrase_concentration'].iloc[0]

# And save
df_environment.to_csv('Salmonella_environment.csv', index=False, sep=',')
