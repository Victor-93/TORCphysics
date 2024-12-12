import sys
#sys.path.append("/users/ph1vvb")
from TORCphysics import params
from TORCphysics import parameter_search as ps
import pandas as pd
import numpy as np
import pickle
from hyperopt import tpe, hp, fmin, Trials

# Description
# ----------------------------------------------------------------------------------------------------------------------
# Same optimization but comparing distributions instead of fitting averages.
# This version _v2 does the same optimization process but it does not vary lacI activity and tetA has a different
# gamma as it induces more supercoils than regular plasmids.

# TODO: We still need to check if the lacI leakage works and does not introduce addiitonal supercoils

# Model & simulation conditions
# ----------------------------------------------------------------------------------------------------------------------
dt = 1.0
final_time = 200 #2000

Ecoli_topoI_C = 17.0
Ecoli_gyrase_C = 44.6
Ecoli_lacI_C = 1.0  # We do not know, at the moment we will assume that lacI behaves the same in Ecoli and in Salmonella
                    # And they both have the same concentration. NOTICE that if concentration=1, I can make the concentration
                    # In Salmonella vary relative to this one, by increasing to 2.0 (twice) or 0.5 (half). Keep that in mind in case
                    # we need it

n_simulations = 12#24 #180 - ask for 61 cores (Each would run 3 simulations per system approximately).
n_batches = 12 # The number of simulations n_simulations will be grouped into n_batches. For each batch, an average
               # will be calculated. In this way, we will be comparing averages with averages.

tests = 4#10 #2000 #1
sigma0 = -0.049 # Initial superhelical density

#lacI_blocking = True  # True if we consider that lacI can block the binding in the minimal PleuWT, False if not.
lacI_blocking = False

#file_out = 'dist_op_TORC_plasmid'
file_out = 'test-dist_op_TORC_plasmid_v2'

if lacI_blocking == True:
    file_out = 'block-'+file_out

ref_file = '../../Experimental_data/reference_distribution.csv'

reference_system = 'Sal_min_WT'
reporter = 'PleuWT'

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

# antitet
oname = 'antitet'
otype = 'site'
presponse = pd.read_csv('../responses/antitet.csv')
binding_model = 'GaussianBinding'
binding_oparams = {
    'k_on': 0.05, 'superhelical_op': -0.06, 'spread': 0.05,  # Parameters for the gaussian binding
    'k_off': 0.05, 'k_closed': 0.02, 'k_open':0.1, 'k_ini':0.3,  # Parameters between stages
    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
variation_dict['antitet'] = variation

# bla
oname = 'bla'
otype = 'site'
presponse = pd.read_csv('../responses/bla.csv')
binding_model = 'GaussianBinding'
binding_oparams = {
    'k_on': 0.05, 'superhelical_op': -0.06, 'spread': 0.05,  # Parameters for the gaussian binding
    'k_off': 0.05, 'k_closed': 0.02, 'k_open':0.1, 'k_ini':0.3,  # Parameters between stages
    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
variation_dict['bla'] = variation

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

# ----------------------------------------------------------------------------------------------------------------------
# Variations
# ----------------------------------------------------------------------------------------------------------------------

# Topoisomerase concentration
# --------------------------------
topoI_C_min = Ecoli_topoI_C/2  # This variation will be applied to Salmonella
topoI_C_max = Ecoli_topoI_C*2
gyrase_C_min = Ecoli_gyrase_C/2
gyrase_C_max = Ecoli_gyrase_C*2

# For genes
# --------------------------------
k_on_min = 0.01
k_on_max = .5
superhelical_op_min=-0.1
superhelical_op_max=0.0
spread_min=0.005
spread_max=0.05
k_closed_min = 0.01
k_closed_max = 0.5
k_open_min = 0.01
k_open_max = 0.5
k_ini_min = 0.01
k_ini_max = 0.5
k_off_min = 0.01
k_off_max = 0.5

# ----------------------------------------------------------------------------------------------------------------------
# PARAMETER SPACE FOR OPTIMIZATION
# ----------------------------------------------------------------------------------------------------------------------
trials = Trials()
space = {
    # topoisomerases
    'topoI_concentration': hp.uniform('topoI_concentration', topoI_C_min, topoI_C_max),
    'gyrase_concentration': hp.uniform('gyrase_concentration', gyrase_C_min, gyrase_C_max),

    # PleuWT
    'PleuWT_k_on': hp.uniform('PleuWT_k_on', k_on_min, k_on_max),
    'PleuWT_spread': hp.uniform('PleuWT_spread', spread_min, spread_max),
    'PleuWT_superhelical_op': hp.uniform('PleuWT_superhelical_op', superhelical_op_min, superhelical_op_max),
    'PleuWT_k_closed': hp.uniform('PleuWT_k_closed', k_closed_min, k_closed_max),
    'PleuWT_k_open': hp.uniform('PleuWT_k_open', k_open_min, k_open_max),
    'PleuWT_k_ini': hp.uniform('PleuWT_k_ini', k_ini_min, k_ini_max),
    'PleuWT_k_off': hp.uniform('PleuWT_k_off', k_off_min, k_off_max),

    # tetA
    'tetA_k_on': hp.uniform('tetA_k_on', k_on_min, k_on_max),
    'tetA_spread': hp.uniform('tetA_spread', spread_min, spread_max),
    'tetA_superhelical_op': hp.uniform('tetA_superhelical_op', superhelical_op_min, superhelical_op_max),
    'tetA_k_closed': hp.uniform('tetA_k_closed', k_closed_min, k_closed_max),
    'tetA_k_open': hp.uniform('tetA_k_open', k_open_min, k_open_max),
    'tetA_k_ini': hp.uniform('tetA_k_ini', k_ini_min, k_ini_max),
    'tetA_k_off': hp.uniform('tetA_k_off', k_off_min, k_off_max),

    # anitet -not at the moment
    'antitet_k_on': hp.uniform('antitet_k_on', k_on_min, k_on_max),
    'antitet_spread': hp.uniform('antitet_spread', spread_min, spread_max),
    'antitet_superhelical_op': hp.uniform('antitet_superhelical_op', superhelical_op_min, superhelical_op_max),
    'antitet_k_closed': hp.uniform('antitet_k_closed', k_closed_min, k_closed_max),
    'antitet_k_open': hp.uniform('antitet_k_open', k_open_min, k_open_max),
    'antitet_k_ini': hp.uniform('antitet_k_ini', k_ini_min, k_ini_max),
    'antitet_k_off': hp.uniform('antitet_k_off', k_off_min, k_off_max),

    # bla - not at the moment
    'bla_k_on': hp.uniform('bla_k_on', k_on_min, k_on_max),
    'bla_spread': hp.uniform('bla_spread', spread_min, spread_max),
    'bla_superhelical_op': hp.uniform('bla_superhelical_op', superhelical_op_min, superhelical_op_max),
    'bla_k_closed': hp.uniform('bla_k_closed', k_closed_min, k_closed_max),
    'bla_k_open': hp.uniform('bla_k_open', k_open_min, k_open_max),
    'bla_k_ini': hp.uniform('bla_k_ini', k_ini_min, k_ini_max),
    'bla_k_off': hp.uniform('bla_k_off', k_off_min, k_off_max)

}

# Circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------
min_circuit_file = '../circuit_min-linear.csv'
complete_circuit_file = '../circuit_complete-linear.csv'

if lacI_blocking:  # If we consider blocking, then we use the displaced lac1 site file which moves the lacO1 a bit
                   # so it can block the promoter
    min_sites_file = 'disp-sites_min-linear_v2.csv'
    complete_sites_file = 'disp-sites_complete-linear_v2.csv'
else:
    min_sites_file = 'sites_min-linear_v2.csv'
    complete_sites_file = 'sites_complete-linear_v2.csv'
enzymes_filename = None
environment_filename = 'environment_v2.csv'
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

    # antitet - site
    antitet_variation = variation_dict['antitet']
    antitet_variation['binding_oparams']['k_on'] = params['antitet_k_on']
    antitet_variation['binding_oparams']['k_off'] = params['antitet_k_off']
    antitet_variation['binding_oparams']['k_open'] = params['antitet_k_open']
    antitet_variation['binding_oparams']['k_closed'] = params['antitet_k_closed']
    antitet_variation['binding_oparams']['k_ini'] = params['antitet_k_ini']
    antitet_variation['binding_oparams']['superhelical_op'] = params['antitet_superhelical_op']
    antitet_variation['binding_oparams']['spread'] = params['antitet_spread']

    # bla - site
    bla_variation = variation_dict['bla']
    bla_variation['binding_oparams']['k_on'] = params['bla_k_on']
    bla_variation['binding_oparams']['k_off'] = params['bla_k_off']
    bla_variation['binding_oparams']['k_open'] = params['bla_k_open']
    bla_variation['binding_oparams']['k_closed'] = params['bla_k_closed']
    bla_variation['binding_oparams']['k_ini'] = params['bla_k_ini']
    bla_variation['binding_oparams']['superhelical_op'] = params['bla_superhelical_op']
    bla_variation['binding_oparams']['spread'] = params['bla_spread']

    # topoI - Environment
    topoI_variation = variation_dict['topoI']
    topoI_variation['concentration'] = params['topoI_concentration']


    # gyrase - Environment
    gyrase_variation = variation_dict['gyrase']
    gyrase_variation['concentration'] = params['gyrase_concentration']

    # Variations independant of calibration (params) - These variations will help us simulate different system
    # conditions

    # ΔtopA - environment
    dtopA_variation = {'name': 'topoI', 'object_type': 'environmental', 'concentration': 0.0}

    # lacI - environment
    dlacI_variation = {'name': 'lacI', 'object_type': 'environmental', 'concentration': 0.0}

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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation ]

    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') & (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'relative'].to_numpy()
    ref_std_value = None
    #ref_std_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') & (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'std'].values[0]

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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation]
    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') & (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'WT'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation, dtopA_variation]

    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation, dtopA_variation]
    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation, dtopA_variation, dlacI_variation]

    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔlacIZYA::FRT\nΔtopA::cat'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation,  dtopA_variation, dlacI_variation]
    ref_value = df.loc[(df['bacterium'] == 'Escherichia coli K12 MG1655') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔlacIZYA::FRT\nΔtopA::cat'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation,
                    topoI_variation, gyrase_variation, dlacI_variation]

    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'WT'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation,
                    topoI_variation, gyrase_variation, dlacI_variation]
    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'WT'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation,
                    gyrase_variation, dlacI_variation, dtopA_variation]

    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].to_numpy()

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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation,
                    gyrase_variation, dlacI_variation, dtopA_variation]
    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔtopA::cat'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation,
                    gyrase_variation,  dtopA_variation]

    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1min mhYFP') & (df['strain'] == 'ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat'), 'relative'].to_numpy()
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
    my_variation = [PleuWT_variation, tetA_variation, antitet_variation, bla_variation,
                    gyrase_variation, dtopA_variation]
    ref_value = df.loc[(df['bacterium'] == 'Salmonella enterica Typhimurium SL1344') &
                       (df['promoter'] == 'PleuWT.1 mhYFP') & (df['strain'] == 'ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat'), 'relative'].to_numpy()
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

#    my_objective, output_dict = ps.KS_calibrate_w_rate(info_list, target_dict, n_simulations,
 #                                                     additional_results=additional_results)
    # my_objective, output_dict = ps.KS_calibrate_w_rate_batches(info_list, target_dict, n_simulations, n_batches,
    #                                                  additional_results=additional_results)
    my_objective, output_dict = ps.KS_calibrate_dist_batches_lossdict(info_list, target_dict, n_simulations, n_batches,
                                                      additional_results=additional_results)


    # Return objective
    if calibrating:
        return my_objective
    else:
        return my_objective, output_dict


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# OUTPUTS
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# Save the current standard output
original_stdout = sys.stdout
# Define the file where you want to save the output
output_file_path = file_out + '.info'

# Open the file in write mode
with open(output_file_path, 'w') as f:
    # Redirect the standard output to the file
    sys.stdout = f

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of TORC Plasmid.")
    print("Launching calibration for dt="+str(dt))
    #print('Ran ' + str(n_simulations) + ' simulations per test. ')
    print('Number of tests = ' + str(tests))

    best = fmin(
        fn=objective_function,  # Objective Function to optimize
        space=space,  # Hyperparameter's Search Space
        algo=tpe.suggest,  # Optimization algorithm (representative TPE)
        max_evals=tests,  # Number of optimization attempts
        trials=trials
    )

    print(" ")
    print("Optimal parameters found from random search: ")
    print(best)

best_df = pd.DataFrame.from_dict([best])
best_df.to_csv(file_out + '.csv', index=False, sep=',')

# Let's save the trials object
# --------------------------------------------------------------------------
with open(file_out+'-trials.pkl', 'wb') as file:
    pickle.dump(trials, file)

# Let's save trials info (params and loses)
# --------------------------------------------------------------------------
params_df = pd.DataFrame(columns=['test', 'loss'])
for n in range(tests):
    tdi = trials.trials[n]  # dictionary with results for test n
    lo = trials.trials[n]['result']['loss']  # loss
    va = trials.trials[n]['misc']['vals']  #values

    # Add a new row using append method
    va['test'] = n
    va['loss'] = lo
    new_row = pd.DataFrame(va)
    params_df = pd.concat([params_df, new_row], ignore_index=True)

params_df.to_csv(file_out + '-values.csv', index=False, sep=',')

# Let's run the function once more with the best params to produce the data so we can then just plot it.
# --------------------------------------------------------------------------
objective, output = objective_function(params=best, calibrating=False)

# Save the dictionary to a file
with open(file_out + '.pkl', 'wb') as file:
    pickle.dump(output, file)
