from TORCphysics import Circuit, Enzyme
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import pandas as pd

# Description: In this script, we will use the Circuit module to produce many simulations of the Pleu500 circuit.
# This one is composed by by the tetA gene transcribed at constant rate which then promotes the transcription of
# the mKalama1 and the Raspberry gene when the bridge formed by the lac proteins is not formed.
# When the bridge is formed, the comunication with mKalama1 is lost

# Initial conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'output'
frames = 10000
series = True
continuation = False
tm = 'continuum'
mm = 'uniform'
dt = 1.0
n_simulations = 1
bridge_time = 5000

for ns in range(n_simulations):
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt, tm, mm)

    my_circuit.name = my_circuit.name + '_' + str(ns)
    my_circuit.log.name = my_circuit.name
    my_circuit.print_general_information()

    Lac1 = Enzyme(e_type=my_circuit.environmental_list[-1].enzyme_type,
                  name=my_circuit.environmental_list[-1].name, site=my_circuit.site_list[5],
                  position=my_circuit.site_list[5].start,
                  size=my_circuit.environmental_list[-1].size, twist=0.0, superhelical=0.0)

    Lac2 = Enzyme(e_type=my_circuit.environmental_list[-1].enzyme_type,
                  name=my_circuit.environmental_list[-1].name, site=my_circuit.site_list[7],
                  position=my_circuit.site_list[7].start + 500,
                  size=my_circuit.environmental_list[-1].size, twist=0.0, superhelical=0.0)

#    for site in my_circuit.site_list:
#        site.k_min = 0.01
    # This is similar to the Run function... but the idea is that we will control when the bridge is formed
    for frame in range(1, frames + 1):
        # print(frame)
        my_circuit.frame = frame
        my_circuit.time = frame * dt
        if my_circuit.series:
            my_circuit.append_sites_to_dict_step1()

        # Apply binding model and get list of new enzymes
        new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt, my_circuit.rng)
        if frame == bridge_time:  # Manually add the lacs
            new_enzyme_list.append(Lac1)
            new_enzyme_list.append(Lac2)
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
        my_circuit.environmental_df = pd.DataFrame.from_dict(my_circuit.environmental_dict_list)
        my_circuit.environmental_df.to_csv(my_circuit.name + '_environment_df.csv', index=False, sep=',')

    # Output the log of events
    my_circuit.log.log_out()

    # Output csvs
    my_circuit.enzyme_list_to_df().to_csv(my_circuit.name + '_enzymes_' + my_circuit.output_prefix + '.csv', index=False, sep=',')
    my_circuit.site_list_to_df().to_csv(my_circuit.name + '_sites_' + my_circuit.output_prefix + '.csv', index=False, sep=',')
    my_circuit.environmental_list_to_df().to_csv(my_circuit.name + '_environment_' + my_circuit.output_prefix + '.csv', index=False, sep=',')
