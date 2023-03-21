from TORCphysics import Circuit, Enzyme
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import pandas as pd


def run_simulation_two_genes(my_circuit):
    RNAP1 = Enzyme(e_type=my_circuit.environmental_list[-1].enzyme_type,
                   name=my_circuit.environmental_list[-1].name, site=my_circuit.site_list[2],
                   position=my_circuit.site_list[2].start,
                   size=my_circuit.environmental_list[-1].size, twist=0.0, superhelical=0.0)

    RNAP2 = Enzyme(e_type=my_circuit.environmental_list[-1].enzyme_type,
                   name=my_circuit.environmental_list[-1].name, site=my_circuit.site_list[3],
                   position=my_circuit.site_list[3].start,
                   size=my_circuit.environmental_list[-1].size, twist=0.0, superhelical=0.0)

    # Let's turn off topoisomerase activity
    my_circuit.environmental_list[0].k_cat = 0.0
    my_circuit.environmental_list[1].k_cat = 0.0

    # This is similar to the Run function... but the idea is that we will control when the bridge is formed
    for frame in range(1, frames + 1):
        my_circuit.frame = frame
        my_circuit.time = frame * dt
        if my_circuit.series:
            my_circuit.append_sites_to_dict_step1()

        # Apply binding model and get list of new enzymes
        if frame == add_time:
            new_enzyme_list = [RNAP1, RNAP2]
        else:
            new_enzyme_list = []
        my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

        # EFFECT
        # --------------------------------------------------------------
        effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt,
                                       my_circuit.topoisomerase_model, my_circuit.mechanical_model)
        my_circuit.apply_effects(effects_list)

        # UNBINDING
        drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list)
        my_circuit.drop_enzymes(drop_list_index)
        my_circuit.add_to_environment(drop_list_enzyme)

        # UPDATE GLOBALS
        my_circuit.update_global_twist()
        my_circuit.update_global_superhelical()

        # Add to series df if the series option was selected (default=True)
        if series:
            my_circuit.append_enzymes_to_dict()
            my_circuit.append_sites_to_dict_step2(new_enzyme_list, drop_list_enzyme)

    # Output the dataframes: (series)
    if series:
        my_circuit.enzymes_df = pd.DataFrame.from_dict(my_circuit.enzymes_dict_list)
        my_circuit.enzymes_df.to_csv(my_circuit.name + '_enzymes_df.csv', index=False, sep=',')
        my_circuit.sites_df = pd.DataFrame.from_dict(my_circuit.sites_dict_list)
        my_circuit.sites_df.to_csv(my_circuit.name + '_sites_df.csv', index=False, sep=',')

    # Output the log of events
    my_circuit.log.log_out()


# Description: Here we will manually add RNAPs to the gene start sites so we can visualize if they behave correctly

# Initial conditions
circuit_filename_linear = 'circuit_linear.csv'
circuit_filename_circular = 'circuit_linear.csv'
sites_filename_tandem = 'sites_tandem.csv'
sites_filename_convergent = 'sites_convergent.csv'
sites_filename_divergent = 'sites_divergent.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'output'
frames = 50
series = True
continuation = False
tm = 'continuum'
mm = 'uniform'
dt = 1.0
add_time = 10  # time in which RNAPs are added

# Let's build the circuits
circuit_lin_tan = Circuit(circuit_filename_linear, sites_filename_tandem, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt, tm, mm)
circuit_lin_dir = Circuit(circuit_filename_linear, sites_filename_divergent, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt, tm, mm)
circuit_lin_con = Circuit(circuit_filename_linear, sites_filename_convergent, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt, tm, mm)

# ---------------------------------------------------------
# Tandem case
# ---------------------------------------------------------
run_simulation_two_genes(circuit_lin_tan)

# TODO: AQUI ME QUEDE, intenta poner la funcion de animacion aqui...
