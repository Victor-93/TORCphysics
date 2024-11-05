from TORCphysics import Circuit, params
import pandas as pd
from TORCphysics import visualization as vs
import numpy as np

# Description
# ----------------------------------------------------------------------------------------------------------------------
# Here, I want to test the model loading and customising for running the random search that will explore which sets
# of parameters or conditions produce the experimental behaviour of the TORC plasmids.

# Model & simulation conditions
# ----------------------------------------------------------------------------------------------------------------------
dt = 1.0
final_time = 500

# ----------------------
# Models to use
# ----------------------
load_list = []

# Sites
# -------------------------------------
# PleuWT
oname = 'PleuWTmin'
otype = 'site'
binding_model = 'SpacerStagesBinding'
presponse = pd.read_csv('../responses/PleuWT.csv')
binding_oparams = {
    'k_on': 0.02, 'superhelical_op': -0.06, 'spacer': 17,  # Parameters for the spacer length binding
    'k_off': 0.1, 'k_closed': 0.02, 'k_open':0.04, 'k_ini':0.1,  # Parameters between stages
    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
load_list.append(variation)

# tetA
oname = 'tetA'
otype = 'site'
presponse = pd.read_csv('../responses/tetA.csv')
binding_model = 'SpacerStagesBinding'
binding_oparams = {
    'k_on': 0.1, 'superhelical_op': -0.06, 'spacer': 17,  # Parameters for the spacer length binding
    'k_off': 0.05, 'k_closed': 0.02, 'k_open':0.1, 'k_ini':0.3,  # Parameters between stages
    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
load_list.append(variation)

# Lacs
oname = 'lac1'
otype = 'site'
binding_model = 'PoissonBinding'
binding_oparams = {'k_on': 0.05}
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
load_list.append(variation)

oname = 'lac2'
variation = {'name': oname, 'object_type': otype,
             'binding_model_name': binding_model, 'binding_oparams': binding_oparams}
load_list.append(variation)

# Environmentals
# -------------------------------------
# topoI
load_list.append({'name': 'topoI', 'object_type': 'environmental', 'concentration': 15.0})

# gyrase
load_list.append({'name': 'gyrase', 'object_type': 'environmental', 'concentration': 40.0})

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
load_list.append(variation)

# Lac - we'll keep the concentration at 1 because we do not know the on/off yet
oname = 'lacI'
otype = 'environmental'
effect_model = 'LacIPoissonBridging'
effect_oparams =  {'k_bridge_on': 0.05, 'k_bridge_off': 0.05, 'leakage': 0.2}
unbinding_model = 'LacISimpleUnBinding'
unbinding_oparams = {'k_off':0.1}
variation = {'name': oname, 'object_type': otype,
             'effect_model_name': effect_model, 'effect_oparams': effect_oparams,
             'unbinding_model_name': unbinding_model, 'unbinding_oparams': unbinding_oparams}
load_list.append(variation)

# Circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------
circuit_filename = '../circuit_min-linear.csv'
sites_filename = 'sites_min-linear.csv'
enzymes_filename = None
environment_filename = 'environment.csv'
output_prefix = 'test'
series = True
continuation = False

# Arange some needed params
initial_time = 0
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Simulation
# ----------------------------------------------------------------------------------------------------------------------
# Load base version of circuit with all sites, environmentals, but without models
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
my_circuit.apply_local_variations(load_list)
my_circuit.print_general_information()
enzymes_df, sites_df, environmental_df = my_circuit.run_return_dfs()

# Animation
# ----------------------------------------------------------------------------------------------------------------------
colors_dict = {'tetA': 'brown', 'PleuWTmin': 'yellow', 'bla': 'purple', 'antitet': 'pink'}
output = 'TORCmin'
out_format = '.mp4'

vs.create_animation_linear(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, output, out_format,
                                  site_type='gene', site_colours=colors_dict)#                           site_type='gene', site_colours=colors_dict)

#my_circuit.run()
