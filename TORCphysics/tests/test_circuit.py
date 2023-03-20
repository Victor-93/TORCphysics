from unittest import TestCase
from TORCphysics import Circuit
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import pandas as pd


class TestCircuit(TestCase):

    # For this test, Circuit should be able to process these non-empty input files (which hopefully are correct).
    # Then, it performs three tests: it counts the number of enzymes, sites and environmentals.
    def test_Circuit_init(self):
        circuit_filename = '../circuit.csv'
        sites_filename = '../sites.csv'
        enzymes_filename = '../enzymes.csv'
        environment_filename = '../environment.csv'
        output_prefix = 'output'
        frames = 5
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        self.assertGreater(my_circuit.get_num_enzymes(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_sites(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_environmentals(), 0, "Empty enzyme list")

    # TODO: This can be several tests in one? 1.- Prove that statistical binding works, 2.- That the effects take place
    #  3.- And that the proteins unbind at the end?
    def test_run(self):
        circuit_filename = '../circuit.csv'
        sites_filename = '../sites.csv'
        enzymes_filename = '../enzymes.csv'
        environment_filename = '../environment.csv'
        output_prefix = 'output'
        frames = 1500
        series = True
        continuation = False
        dt = 0.5
        tm = 'continuum'
        mm = 'uniform'
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        #        for site in my_circuit.site_list:  # I'll increase the rates
        #             site.k_min = site.k_min * 100
        # site.k_min = site.k_min * 5
        # Let's make the rates
        my_circuit.run()
        print(0)

    def test_run2(self):
        circuit_filename = '../circuit.csv'
        sites_filename = '../sites.csv'
        enzymes_filename = '../enzymes.csv'
        environment_filename = '../environment.csv'
        output_prefix = 'output'
        frames = 1500
        series = True
        continuation = False
        dt = 0.5
        tm = 'continuum'
        mm = 'uniform'
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.environmental_list[0].concentration = 0
        my_circuit.environmental_list[1].concentration = 0
        print(my_circuit.environmental_list)
        for site in my_circuit.site_list:  # I'll increase the rates
            # site.k_min = site.k_min * 100
            site.k_min = site.k_min * 0
        print(my_circuit.site_list[2].name)
        my_circuit.site_list[2].k_min = 0.01
        # Let's make the rates
        my_circuit.run()
        print(0)

    # TODO: test fake ends - circular and linear. No enzymes, 1 enzyme, 2 enzymes
    # Tests the positions of the fake enzymes
    def test_fake_ends(self):
        # CIRCULAR
        # -----------------------------------------------------------------------------------------
        # Empty case
        circuit_filename = 'test_inputs/test_circuit_circular.csv'
        sites_filename = 'test_inputs/sites_1_gene.csv'
        enzymes_filename = 'test_inputs/empty_enzymes.csv'
        environment_filename = 'test_inputs/RNAP_environment.csv'
        output_prefix = 'output'
        frames = 5
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        self.assertEqual(my_circuit.enzyme_list[0].position, 0, "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit.enzyme_list[-1].position, my_circuit.size + 1,
                         "Wrong position of right fake enzyme")

        # One Enzyme bound
        enzymes1_filename = 'test_inputs/enzymes_1.csv'
        my_circuit_1 = Circuit(circuit_filename, sites_filename, enzymes1_filename, environment_filename,
                               output_prefix, frames, series, continuation, dt, tm, mm)
        enzyme1 = my_circuit_1.enzyme_list[1]
#        self.assertEqual(my_circuit_1.enzyme_list[0].position, 1 + enzyme1.position + enzyme1.size - my_circuit_1.size,
#                         "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit_1.enzyme_list[0].position,  enzyme1.position + enzyme1.size - my_circuit_1.size,
                         "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit_1.enzyme_list[-1].position, enzyme1.position + my_circuit_1.size,
                         "Wrong position of right fake enzyme")

        # Two Enzymes bound
        enzymes2_filename = 'test_inputs/enzymes_2.csv'
        my_circuit_2 = Circuit(circuit_filename, sites_filename, enzymes2_filename, environment_filename,
                               output_prefix, frames, series, continuation, dt, tm, mm)
        enzyme1 = my_circuit_2.enzyme_list[1]
        enzyme2 = my_circuit_2.enzyme_list[2]
#        self.assertEqual(my_circuit_2.enzyme_list[0].position, 1 + enzyme2.position + enzyme2.size - my_circuit_2.size,
#                         "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit_2.enzyme_list[0].position, enzyme2.position + enzyme2.size - my_circuit_2.size,
                         "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit_2.enzyme_list[-1].position, enzyme1.position + my_circuit_2.size,
                         "Wrong position of right fake enzyme")
        # LINEAR
        # -----------------------------------------------------------------------------------------
        circuit_filename_linear = 'test_inputs/test_circuit_linear.csv'
        # No enzymes
        my_circuit_linear = Circuit(circuit_filename_linear, sites_filename, enzymes_filename, environment_filename,
                                    output_prefix, frames, series, continuation, dt, tm, mm)
        self.assertEqual(my_circuit_linear.enzyme_list[0].position, 0, "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit_linear.enzyme_list[-1].position, my_circuit_linear.size + 1,
                         "Wrong position of right fake enzyme")
        # 1 enzyme
        my_circuit_linear1 = Circuit(circuit_filename_linear, sites_filename, enzymes1_filename, environment_filename,
                                     output_prefix, frames, series, continuation, dt, tm, mm)
        self.assertEqual(my_circuit_linear1.enzyme_list[0].position, 0, "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit_linear1.enzyme_list[-1].position, my_circuit_linear1.size + 1,
                         "Wrong position of right fake enzyme")
        # 2 enzymes
        my_circuit_linear2 = Circuit(circuit_filename_linear, sites_filename, enzymes2_filename, environment_filename,
                                     output_prefix, frames, series, continuation, dt, tm, mm)
        self.assertEqual(my_circuit_linear2.enzyme_list[0].position, 0, "Wrong position of left fake enzyme")
        self.assertEqual(my_circuit_linear2.enzyme_list[-1].position, my_circuit_linear2.size + 1,
                         "Wrong position of right fake enzyme")

    # TODO: test bind/unbind no topo - one gene, a RNAP binds, then unbinds and twist is conserved. - test with 2 genes
    def test_bind_unbind_1_RNAP_1_gene(self):
        circuit_filename = 'test_inputs/test_circuit_circular.csv'
        sites_filename = 'test_inputs/sites_1_gene.csv'
        enzymes_filename = 'test_inputs/empty_enzymes.csv'
        environment_filename = 'test_inputs/RNAP_environment.csv'
        output_prefix = 'output'
        frames = 500
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.site_list[2].k_min = 0.0
        s0 = my_circuit.superhelical

        # This is similar to the Run function... but the idea is that we will control the rate
        for frame in range(1, frames + 1):
            my_circuit.frame = frame
            my_circuit.time = frame * dt
            if my_circuit.series:
                my_circuit.append_sites_to_dict_step1()

            if frame == 100:
                my_circuit.site_list[2].k_min = .1

            # Apply binding model and get list of new enzymes
            new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt)

            if len(new_enzyme_list) > 0:
                my_circuit.site_list[2].k_min = 0.0
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
        sf = my_circuit.superhelical  # final superhelical
        self.assertEqual(abs(s0-sf), 0, "Superhelical changed")

    # Bind unbind with 1 NAP on the right
    def test_bind_unbind_1_RNAP_1_gene_1_NAP_right(self):
        circuit_filename = 'test_inputs/test_circuit_circular.csv'
        sites_filename = 'test_inputs/sites_1_gene.csv'
        enzymes_filename = 'test_inputs/enzymes_1.csv'
        environment_filename = 'test_inputs/RNAP_environment.csv'
        output_prefix = 'output'
        frames = 500
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.site_list[2].k_min = 0.0
        s0 = my_circuit.superhelical
        t0 = my_circuit.twist
        err = .00000000000000001

        # This is similar to the Run function... but the idea is that we will control the rate
        for frame in range(1, frames + 1):
            my_circuit.frame = frame
            my_circuit.time = frame * dt
            if my_circuit.series:
                my_circuit.append_sites_to_dict_step1()

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error initial, frame', frame)

            if frame == 10:
                my_circuit.site_list[2].k_min = .01

            if frame == 450:
                my_circuit.site_list[2].k_min = .0

            # Apply binding model and get list of new enzymes
            new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt)

#            if len(new_enzyme_list) > 0:
#                my_circuit.site_list[2].k_min = 0.0
            my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error new_enzymes frame', frame)
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

            if abs(t0-my_circuit.twist) > err:
                print('error drop_enzymes, frame', frame)


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
        sf = my_circuit.superhelical  # final superhelical
        self.assertLessEqual(abs(s0-sf), err, "Superhelical changed")

    # Bind unbind with 1 NAP on the left
    def test_bind_unbind_1_RNAP_1_gene_1_NAP_left(self):
        circuit_filename = 'test_inputs/test_circuit_circular.csv'
        sites_filename = 'test_inputs/sites_1_gene_4500.csv'
        enzymes_filename = 'test_inputs/enzymes_1.csv'
        environment_filename = 'test_inputs/RNAP_environment.csv'
        output_prefix = 'output'
        frames = 500
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.site_list[2].k_min = 0.0
        s0 = my_circuit.superhelical
        t0 = my_circuit.twist
        err = .01

        # This is similar to the Run function... but the idea is that we will control the rate
        for frame in range(1, frames + 1):
            my_circuit.frame = frame
            my_circuit.time = frame * dt
            if my_circuit.series:
                my_circuit.append_sites_to_dict_step1()

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error initial, frame', frame)

            if frame == 100:
                my_circuit.site_list[2].k_min = .1

            # Apply binding model and get list of new enzymes
            new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt)

            if len(new_enzyme_list) > 0:
                my_circuit.site_list[2].k_min = 0.0
            my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error new_enzymes frame', frame)
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

            if abs(t0-my_circuit.twist) > err:
                print('error drop_enzymes, frame', frame)


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
        sf = my_circuit.superhelical  # final superhelical
        self.assertLessEqual(abs(s0-sf), err, "Superhelical changed")

    # Bind unbind with 1 NAP on the left and right
    def test_bind_unbind_1_RNAP_1_gene_1_NAP_left_right(self):
        circuit_filename = 'test_inputs/test_circuit_circular.csv'
        sites_filename = 'test_inputs/sites_1_gene_4500.csv'
        enzymes_filename = 'test_inputs/enzymes_2.csv'
        environment_filename = 'test_inputs/RNAP_environment.csv'
        output_prefix = 'output'
        frames = 500
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.site_list[2].k_min = 0.0
        s0 = my_circuit.superhelical
        t0 = my_circuit.twist
        err = .1

        # This is similar to the Run function... but the idea is that we will control the rate
        for frame in range(1, frames + 1):
            my_circuit.frame = frame
            my_circuit.time = frame * dt
            if my_circuit.series:
                my_circuit.append_sites_to_dict_step1()

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error initial, frame', frame)

            if frame == 100:
                my_circuit.site_list[2].k_min = .1

            if frame == 400:
                my_circuit.site_list[2].k_min = 0.0

            # Apply binding model and get list of new enzymes
            new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt)

            if len(new_enzyme_list) > 0:
                my_circuit.site_list[2].k_min = 0.0
            my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error new_enzymes frame', frame)
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

            if abs(t0-my_circuit.twist) > err:
                print('error drop_enzymes, frame', frame)


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
        sf = my_circuit.superhelical  # final superhelical
        self.assertLessEqual(abs(s0-sf), err, "Superhelical changed")


    def test_bind_unbind_many_RNAPs_1_gene(self):
        circuit_filename = 'test_inputs/test_circuit_circular.csv'
        sites_filename = 'test_inputs/sites_1_gene.csv'
        enzymes_filename = 'test_inputs/empty_enzymes.csv'
        environment_filename = 'test_inputs/RNAP_environment.csv'
        output_prefix = 'output'
        frames = 1000
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.site_list[2].k_min = 0.0
        s0 = my_circuit.superhelical
        t0 = my_circuit.twist
        err = 0.1
        # This is similar to the Run function... but the idea is that we will control the rate
        for frame in range(1, frames + 1):
            #print(frame)
            my_circuit.frame = frame
            my_circuit.time = frame * dt
            if my_circuit.series:
                my_circuit.append_sites_to_dict_step1()

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error initial, frame', frame)

            if frame == 100:
                my_circuit.site_list[2].k_min = .05
            elif frame == 900:
                my_circuit.site_list[2].k_min = 0.0

            # Apply binding model and get list of new enzymes
            new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt)

            my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            t2 = my_circuit.twist+2
            ss = my_circuit.enzyme_list[0].twist
            s1 = sum(enzyme.twist for enzyme in my_circuit.enzyme_list)
            my_circuit.update_global_twist()
            t3 = my_circuit.twist+3
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('new enzyme, frame', frame)

            # EFFECT
            # --------------------------------------------------------------
            effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt,
                                           my_circuit.topoisomerase_model, my_circuit.mechanical_model)
            my_circuit.apply_effects(effects_list)

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('effects, frame', frame)

            # UNBINDING
            drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list)
            my_circuit.drop_enzymes(drop_list_index)

            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('error initial, frame', frame)
            my_circuit.add_to_environment(drop_list_enzyme)

            # UPDATE GLOBALS
            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            if abs(t0-my_circuit.twist) > err:
                print('Adding to the environment, frame', frame)

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
        sf = my_circuit.superhelical  # final superhelical
        self.assertLessEqual(abs(s0-sf), err, "Superhelical changed")

    # TODO: test empty circuit -

#    def test_add_1_enzyme(self):
