import sys
from TORCphysics import params
from TORCphysics import parameter_search as ps
import pandas as pd
import numpy as np
import pickle

# Description
# ----------------------------------------------------------------------------------------------------------------------
# Following the optimization proces, we want to reproduce the experiment and get some data.
# This script produces a pickle file with data called file_out .pkl (check the file_out variable  in the input).
# This pickle file is a list of dictionaries with the results, and each entry represents a different system e.g., Ecoli WT.
# For each experiment, it provides system information such as name, description, bacterium, promoter and strain.
# It also includes the 'reference' which is the relative expression rate from experiments,
# the 'relative_rate' which is the relative expression rate calculated from simulations, the objective component
# (difference between relative expressions with experimental data), the number of transcripts of reporter gene,
# 'global_superhelical' density, 'local_superhelical' density per site, and the 'production_rate' per site as well.

# INPUT: Model & simulation conditions
# ----------------------------------------------------------------------------------------------------------------------
dt = 1.0
final_time = 10000

Ecoli_topoI_C = 17.0
Ecoli_gyrase_C = 44.6
Ecoli_lacI_C = 1.0  # We do not know, at the moment we will assume that lacI behaves the same in Ecoli and in Salmonella
                    # And they both have the same concentration. NOTICE that if concentration=1, I can make the concentration
                    # In Salmonella vary relative to this one, by increasing to 2.0 (twice) or 0.5 (half). Keep that in mind in case
                    # we need it

n_simulations = 12#50
sigma0 = -0.049 # Initial superhelical density

lacI_blocking = True

file_out='rep-TORC_plasmid'

if lacI_blocking == True:
    file_out = 'block-'+file_out

ref_file = '../../Experimental_data/reference.csv'

reference_system = 'Sal_min_WT'
reporter = 'PleuWT'

param_file = 'calibration_TORC_plasmid.csv'   # This one has the parametrisation from calibration

# For Salmonella we do not know! So we will vary it

# ----------------------------------------------------------------------------------------------------------------------
# Models to use
# ----------------------------------------------------------------------------------------------------------------------
load_list = []
variation_dict = {}

# Sites
# -------------------------------------
# PleuWT
oname = 'PleuWT'
otype = 'site'
binding_model = 'GaussianBinding'
presponse = pd.read_csv('../responses/PleuWT.csv')
binding_oparams = {
    'k_on': 0.05, 'superhelical_op': -0.06, 'spread': 0.05,  # Parameters for the gaussian binding
    'k_off': 0.1, 'k_closed': 0.02, 'k_open':0.04, 'k_ini':0.1,  # Parameters between stages
    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
variation_dict['PleuWT'] = variation

#load_list.append(variation)

# tetA
oname = 'tetA'
otype = 'site'
presponse = pd.read_csv('../responses/tetA.csv')
binding_model = 'GaussianBinding'
binding_oparams = {
    'k_on': 0.05, 'superhelical_op': -0.06, 'spread': 0.05,  # Parameters for the gaussian binding
    'k_off': 0.05, 'k_closed': 0.02, 'k_open':0.1, 'k_ini':0.3,  # Parameters between stages
    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
variation_dict['tetA'] = variation

#load_list.append(variation)

# Lacs
oname = 'lac1'
otype = 'site'
binding_model = 'PoissonBinding'
binding_oparams = {'k_on': 0.05}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
#load_list.append(variation)
variation_dict['lac1'] = variation

oname = 'lac2'
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}

variation_dict['lac2'] = variation

#load_list.append(variation)

# Environmentals
# -------------------------------------
# topoI
variation = {'name': 'topoI', 'object_type': 'environmental', 'concentration': 15.0}
#load_list.append({'name': 'topoI', 'object_type': 'environmental', 'concentration': 15.0})
#load_list.append(variation)
variation_dict['topoI'] = variation

# gyrase
variation = {'name': 'gyrase', 'object_type': 'environmental', 'concentration': 40.0}
#load_list.append({'name': 'gyrase', 'object_type': 'environmental', 'concentration': 40.0})
variation_dict['gyrase'] = variation

# RNAP
oname = 'RNAP'
otype = 'environmental'
effect_model = 'RNAPStagesStallv2'
effect_oparams =  {'velocity': params.v0, 'kappa': params.RNAP_kappa, 'stall_torque': params.stall_torque,
                   'gamma':0.05}
unbinding_model = 'RNAPStagesSimpleUnbindingv2'
variation = {'name': oname, 'object_type': otype,
             'effect_model_name': effect_model, 'effect_oparams': effect_oparams,
             'unbinding_model_name': unbinding_model, 'unbinding_oparams': {}}
variation_dict['RNAP'] = variation
#load_list.append(variation)

# Lac - we'll keep the concentration at 1 because we do not know the on/off yet
oname = 'lacI'
otype = 'environmental'
#binding_model = 'PoissonBinding'
#binding_oparams = {'k_on': 0.2}
effect_model = 'LacIPoissonBridging'
effect_oparams =  {'k_bridge_on': 0.05, 'k_bridge_off': 0.05, 'leakage': 0.2}
unbinding_model = 'LacISimpleUnBinding'
unbinding_oparams = {'k_off':0.1}
variation = {'name': oname, 'object_type': otype,
 #            'binding_model_name': binding_model, 'binding_oparams': binding_oparams,
             'effect_model_name': effect_model, 'effect_oparams': effect_oparams,
             'unbinding_model_name': unbinding_model, 'unbinding_oparams': unbinding_oparams}
variation_dict['lacI'] = variation
#load_list.append(variation)

# ----------------------------------------------------------------------------------------------------------------------
# OPTIMIZATION params - Let's read the csv file
# ----------------------------------------------------------------------------------------------------------------------
calibration_params = pd.read_csv(param_file)

# Convert the DataFrame to a dictionary (list of dictionaries)
calibration_dict = calibration_params.to_dict(orient='records')[0]  # Get the first (and only) dictionary

# Now calibration_dict will have the key 'kon' instead of 'k_on'
print('This is the calibration params: ', calibration_dict)

# Circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------
min_circuit_file = '../circuit_min-linear.csv'
complete_circuit_file = '../circuit_complete-linear.csv'
if lacI_blocking:  # If we consider blocking, then we use the displaced lac1 site file which moves the lacO1 a bit
                   # so it can block the promoter
    min_sites_file = 'disp-sites_min-linear.csv'
    complete_sites_file = 'disp-sites_complete-linear.csv'
else:
    min_sites_file = 'sites_min-linear.csv'
    complete_sites_file = 'sites_complete-linear.csv'
enzymes_filename = None
environment_filename = 'environment.csv'
output_prefix = 'TORC_plasmid'
series = True
continuation = False

# Arange some needed params
initial_time = 0
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Reference
# ----------------------------------------------------------------------------------------------------------------------
# Let's build our reference or info file
reference_pd = pd.read_csv(ref_file)


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Optimization function
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function and can return the outputs
# if the calibration option is not set (reproduction)

def objective_function(params, calibrating=True):

    # ------------------------------------------------------------------------------------------------------------------
    # Apply variations (params) - Let's do this first before building our systems
    # ------------------------------------------------------------------------------------------------------------------
    # PleuWT - site
    PleuWT_variation = variation_dict['PleuWT']
    PleuWT_variation['binding_oparams']['k_on'] = params['PleuWT_k_on']
    PleuWT_variation['binding_oparams']['k_off'] = params['PleuWT_k_off']
    PleuWT_variation['binding_oparams']['k_open'] = params['PleuWT_k_open']
    PleuWT_variation['binding_oparams']['k_closed'] = params['PleuWT_k_closed']
    PleuWT_variation['binding_oparams']['k_ini'] = params['PleuWT_k_ini']
    PleuWT_variation['binding_oparams']['superhelical_op'] = params['PleuWT_superhelical_op']
    PleuWT_variation['binding_oparams']['spread'] = params['PleuWT_spread']

    # tetA - site
    tetA_variation = variation_dict['tetA']
    tetA_variation['binding_oparams']['k_on'] = params['tetA_k_on']
    tetA_variation['binding_oparams']['k_off'] = params['tetA_k_off']
    tetA_variation['binding_oparams']['k_open'] = params['tetA_k_open']
    tetA_variation['binding_oparams']['k_closed'] = params['tetA_k_closed']
    tetA_variation['binding_oparams']['k_ini'] = params['tetA_k_ini']
    tetA_variation['binding_oparams']['superhelical_op'] = params['tetA_superhelical_op']
    tetA_variation['binding_oparams']['spread'] = params['tetA_spread']

    # lac1/lac2 - site
    lac1_variation = variation_dict['lac1']
    lac1_variation['binding_oparams']['k_on'] = params['lacI_k_on']
    lac2_variation = variation_dict['lac2']
    lac2_variation['binding_oparams']['k_on'] = params['lacI_k_on']

    # topoI - Environment
    topoI_variation = variation_dict['topoI']
    topoI_variation['concentration'] = params['topoI_concentration']


    # gyrase - Environment
    gyrase_variation = variation_dict['gyrase']
    gyrase_variation['concentration'] = params['gyrase_concentration']

    # lacI - Environment
    lacI_variation = variation_dict['lacI']
    lacI_variation['effect_oparams']['k_bridge_on'] = params['bridge_on']
    lacI_variation['effect_oparams']['k_bridge_off'] = params['bridge_off']
    lacI_variation['effect_oparams']['leakage'] = params['leakage']
    lacI_variation['unbinding_oparams']['k_off'] = params['lacI_k_off']

    # Variations independant of calibration (params) - These variations will help us simulate different system
    # conditions

    # ΔtopA - environment
    dtopA_variation = {'name': 'topoI', 'object_type': 'environmental', 'concentration': 0.0}

    # lacI - environment
    dlacI_variation = {'name': 'lacI', 'object_type': 'environmental', 'concentration': 0.0}

    # RNAP - Environment
    RNAP_variation = variation_dict['RNAP']


    # ------------------------------------------------------------------------------------------------------------------
    # Build systems - Sort out input information to run each of the systems with their corresponding variations
    # ------------------------------------------------------------------------------------------------------------------

    # We have 6 systems with different conditions and each with two promoters (minimal and complete).
    # So, in total, we have 2*6=12 systems:
    # 1.- EColi WT  (topoI, gyrase and lacI)
    # 2.- EColi ΔtopA  (gyrase and lacI)
    # 3.- Ecoli ΔtopA-ΔlacI (gyrase)
    # 4.- Sal WT (topoI and gyrase)
    # 5.- Sal ΔtopA (gyrase)
    # 6.- Sal ΔtopA-lacI (gyrase and lacI)

    # This variable will contain all the information needed to simulate the systems and will store the calculated
    # rates and relative rates. Each entry will be a dictionary that has all the relevant information for the system
    info_list = []

    df = reference_pd  # This is the relative expression rates from experimental data. I just rename it to make
                       # operations shorter so keep this in mind!

    # =================================================
    # =================================================
    # ECOLI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # =================================================
    # =================================================
    # Probably there is a more compact way to define these systems, but I'll make it more explicitly anyway

    # 1.1.- E. Coli, Minimum, WT (topoI, gyrase and lacI):
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': min_circuit_file, 'sites_filename': min_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      lacI_variation, RNAP_variation]

    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') & (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') & (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'std'].values[0]

    info_dict = {'name':'EColi_min_WT', 'description': ' E. Coli, Minimum promoter, WT background',
                 'bacterium': 'EColi', 'promoter': 'min', 'strain': 'WT',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 1.2.- E. Coli, Complete, WT:
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': complete_circuit_file, 'sites_filename': complete_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      lacI_variation, RNAP_variation]
    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') & (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'WT'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') & (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'WT'), 'std'].values[0]
    info_dict = {'name': 'EColi_full_WT', 'description': ' E. Coli, Complete promoter, WT background',
                 'bacterium': 'EColi', 'promoter': 'full', 'strain': 'WT',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 2.1.- E. Coli, Minimum, ΔtopA (gyrase and lacI):
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': min_circuit_file, 'sites_filename': min_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      lacI_variation, RNAP_variation, dtopA_variation]

    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'std'].values[0]
    info_dict = {'name': 'EColi_min_dtopA', 'description': ' E. Coli, Minimum promoter, ΔtopA background',
                 'bacterium': 'EColi', 'promoter': 'min', 'strain': 'ΔtopA',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 2.2.- E. Coli, Complete, ΔtopA:
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': complete_circuit_file, 'sites_filename': complete_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      lacI_variation, RNAP_variation, dtopA_variation]
    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'std'].values[0]
    info_dict = {'name': 'EColi_full_dtopA', 'description': ' E. Coli, Complete promoter, ΔtopA background',
                 'bacterium': 'EColi', 'promoter': 'full', 'strain': 'ΔtopA',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 3.1.- E. Coli, Minimum, ΔtopA-ΔlacI (gyrase):
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': min_circuit_file, 'sites_filename': min_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, dtopA_variation, dlacI_variation]

    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔlacIZYA::FRT\nΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔlacIZYA::FRT\nΔtopA::cat'), 'std'].values[0]
    info_dict = {'name':'EColi_min_dtopA-dlacI', 'description': ' E. Coli, Minimum promoter, ΔtopA-ΔlacI background',
                 'bacterium': 'EColi', 'promoter': 'min', 'strain': 'ΔtopA-ΔlacI',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 3.2.- E. Coli, Complete, ΔtopA-ΔlacI:
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': complete_circuit_file, 'sites_filename': complete_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, dtopA_variation, dlacI_variation]
    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔlacIZYA::FRT\nΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (
                                   df['strain'] == 'ΔlacIZYA::FRT\nΔtopA::cat'), 'std'].values[0]
    info_dict = {'name': 'EColi_full_dtopA-dlacI','description': ' E. Coli, Complete promoter, ΔtopA-ΔlacI background',
                 'bacterium': 'EColi', 'promoter': 'full', 'strain': 'ΔtopA-ΔlacI',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # =================================================
    # =================================================
    # SALMONELLA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # =================================================
    # =================================================

    # 4.1.- Salmonella, Minimum, WT (topoI, gyrase):
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': min_circuit_file, 'sites_filename': min_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, topoI_variation, gyrase_variation, dlacI_variation]

    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'std'].values[0]
    info_dict = {'name':'Sal_min_WT','description': 'Salmonella, Minimum promoter, WT background',
                 'bacterium': 'Salmonella', 'promoter': 'min', 'strain': 'WT',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 4.2.- Salmonella, Complete, WT:
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': complete_circuit_file, 'sites_filename': complete_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, topoI_variation, gyrase_variation, dlacI_variation]
    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'WT'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'WT'), 'std'].values[0]
    info_dict = {'name': 'Sal_full_WT', 'description': 'Salmonella, Complete promoter, WT background',
                 'bacterium': 'Salmonella', 'promoter': 'full', 'strain': 'WT',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 5.1.- Salmonella, Minimum, ΔtopA (gyrase):
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': min_circuit_file, 'sites_filename': min_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, gyrase_variation, dlacI_variation, dtopA_variation]

    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'std'].values[0]

    info_dict = {'name':'Sal_min_dtopA', 'description': 'Salmonella, Minimum promoter, ΔtopA background',
                 'bacterium': 'Salmonella', 'promoter': 'min', 'strain': 'ΔtopA',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 5.2.- Salmonella, Complete, ΔtopA:
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': complete_circuit_file, 'sites_filename': complete_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, gyrase_variation, dlacI_variation, dtopA_variation]
    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'std'].values[0]
    info_dict = {'name':'Sal_full_dtopA', 'description': 'Salmonella, Complete promoter, ΔtopA background',
                 'bacterium': 'Salmonella', 'promoter': 'full', 'strain': 'ΔtopA',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 6.- Sal ΔtopA-lacI (gyrase and lacI)
    # 6.1.- Salmonella, Minimum, ΔtopA-lacI (gyrase and lacI):
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': min_circuit_file, 'sites_filename': min_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, gyrase_variation, lacI_variation, dtopA_variation]

    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat'), 'std'].values[0]
    info_dict = {'name':'Sal_min_dtopA-lacI','description': 'Salmonella, Minimum promoter, ΔtopA-lacI background',
                 'bacterium': 'Salmonella', 'promoter': 'min', 'strain': 'ΔtopA-lacI',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # 6.2.- Salmonella, Complete, WT:
    # ------------------------------------------------------------------------------------------------------------------
    my_global_dict = {'circuit_filename': complete_circuit_file, 'sites_filename': complete_sites_file,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma0,
                       'DNA_concentration': 0.0}
    my_variation = [PleuWT_variation, tetA_variation, lac1_variation, lac2_variation,
                      RNAP_variation, gyrase_variation, lacI_variation, dtopA_variation]
    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat'), 'relative'].values[0]
    ref_std_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (
                                   df['strain'] == 'ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat'), 'std'].values[0]
    info_dict = {'name':'Sal_full_dtopA-lacI', 'description': 'Salmonella, Complete promoter, ΔtopA-lacI background',
                 'bacterium': 'Salmonella', 'promoter': 'full', 'strain': 'ΔtopA-lacI',
                 'reference': ref_value, 'reference_std': ref_std_value, 'global_conditions': my_global_dict,
                 'variations': my_variation, 'transcripts': None, 'relative_rate': None}
    info_list.append(info_dict)

    # Note that it can be expanded to have multiple reporters (maybe if you treat it as a list)
    target_dict = {'reference_system': reference_system, 'reporter': reporter}

    # Finally, run objective function.
    # ------------------------------------------
    if calibrating:
        additional_results = False
    else:
        additional_results = True

    my_objective, output_dict = ps.calibrate_w_rate(info_list, target_dict, n_simulations,
                                                      additional_results=additional_results)

    # Return objective
    if calibrating:
        return my_objective
    else:
        return my_objective, output_dict


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# REPRODUCE
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Let's run the function once more with the best params to produce the data so we can then just plot it.
objective, output = objective_function(params=calibration_dict, calibrating=False)

# Save the dictionary to a file
with open(file_out + '.pkl', 'wb') as file:
    pickle.dump(output, file)
