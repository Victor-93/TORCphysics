from TORCphysics import Enzyme
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import pandas as pd


class AddEnzyme:
    # I thought it'll be easier to describe the effects as an object.
    # Because the current effects are taken place at the current enzyme i=index, I use the index to locate the enzyme
    # in the enzyme_list which is applying the effect.
    # These effects act locally, so they can modify the enzyme's position, and the twist at the neighbouring domains.
    def __init__(self, frame, e_type, name, site, position, size):
        # I'll save the input filenames just in case#
        self.frame = frame  # the frame to add the enzyme
        # Define enzyme. Do you think we would need the twist and superhelical? Probably no.
        self.enzyme = Enzyme(e_type=e_type, name=name, site=site, position=position, size=size,
                             twist=0.0, superhelical=0.0, k_cat=0.0, k_off=0.0)


# Run normal simulation but adding custom enzymes.
# custom enzymes is a list of CustomEnzyme objects...
def run_simulation(my_circuit, custom_enzymes):
    custom_enzymes.sort(key=lambda x: x.frame)  # Sort custom enzymes by time

    for frame in range(1, my_circuit.frames + 1):
        custom_enzymes_frame = [enzyme for enzyme in custom_enzymes if enzyme.frame == frame]

        my_circuit.frame += 1
        my_circuit.time = frame * my_circuit.dt
        if my_circuit.series:
            my_circuit.append_sites_to_dict_step1()

        # If we are modeling the topoisomerase binding as a stochastic process, we need to define the sites.
        # if my_circuit.topoisomerase_model == 'stochastic':
        #    my_circuit.define_topoisomerase_binding_sites()
        # BINDING
        # --------------------------------------------------------------
        # Apply binding model and get list of new enzymes
        new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                           my_circuit.rng)

        # ADD CUSTOM ENZYMES
        # --------------------------------------------------------------
        if custom_enzymes_frame:  # If list is not empty
            for custom_enzyme in custom_enzymes_frame:
                new_enzyme_list.append(custom_enzyme.enzyme)

        # These new enzymes are lacking twist and superhelical, we need to fix them and actually add them
        # to the enzyme_list
        # But before, add the binding events to the log  (it's easier to do it first)
        #            my_circuit.add_binding_events_to_log(new_enzyme_list)
        my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

        # EFFECT
        # --------------------------------------------------------------
        effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                       my_circuit.topoisomerase_model, my_circuit.mechanical_model)
        my_circuit.apply_effects(effects_list)

        # UNBINDING
        # --------------------------------------------------------------
        drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list, my_circuit.dt, my_circuit.rng)
        my_circuit.drop_enzymes(drop_list_index)
        my_circuit.add_to_environment(drop_list_enzyme)
        #            my_circuit.add_unbinding_events_to_log(drop_list)

        # UPDATE GLOBALS
        # --------------------------------------------------------------
        my_circuit.update_global_twist()
        my_circuit.update_global_superhelical()

        # my_circuit.log.log_out()

        # Add to series df if the series option was selected (default=True)
        # --------------------------------------------------------------
        if my_circuit.series:
            my_circuit.append_enzymes_to_dict()
            my_circuit.append_sites_to_dict_step2(new_enzyme_list, drop_list_enzyme)
            my_circuit.append_environmental_to_dict()

    # Output the dataframes: (series)
    if my_circuit.series:
        my_circuit.enzymes_df = pd.DataFrame.from_dict(my_circuit.enzymes_dict_list)
        my_circuit.enzymes_df.to_csv(my_circuit.name + '_' + my_circuit.output_prefix + '_enzymes_df.csv', index=False,
                                     sep=',')
        my_circuit.sites_df = pd.DataFrame.from_dict(my_circuit.sites_dict_list)
        my_circuit.sites_df.to_csv(my_circuit.name + '_' + my_circuit.output_prefix + '_sites_df.csv', index=False,
                                   sep=',')
        my_circuit.environmental_df = pd.DataFrame.from_dict(my_circuit.environmental_dict_list)
        my_circuit.environmental_df.to_csv(my_circuit.name + '_' + my_circuit.output_prefix + '_environment_df.csv',
                                           index=False, sep=',')

    # Output the log of events
    my_circuit.log.final_twist = my_circuit.twist
    my_circuit.log.final_superhelical = my_circuit.superhelical
    my_circuit.log.log_out()

    # Output csvs
    my_circuit.enzyme_list_to_df().to_csv(my_circuit.name + '_' + my_circuit.output_prefix + '_enzymes' + '.csv',
                                          index=False, sep=',')
    my_circuit.site_list_to_df().to_csv(my_circuit.name + '_' + my_circuit.output_prefix + '_sites' + '.csv',
                                        index=False, sep=',')
    my_circuit.environmental_list_to_df().to_csv(
        my_circuit.name + '_' + my_circuit.output_prefix + '_environment' + '.csv',
        index=False, sep=',')
