import copy
from unittest import TestCase
from TORCphysics import Circuit, Enzyme
from TORCphysics import models_workflow as mw
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import pandas as pd

circuit_circular_filename = 'circuit_circular.csv'
circuit_linear_filename = 'circuit_linear.csv'
sites_filename = 'sites.csv'  # Has one lacOP binding site
enzymes_filename = '../test_inputs/empty_enzymes.csv'
environment_filename = 'environment.csv'  # Has specific binding enzymes and bare DNA binding. All with concentration=0
output_prefix = ''
frames = 1000
series = True
continuation = False
dt = 1
error = 0.0001


# TODO: Maybe my code wasn't wrong, because I wrote that when molecules bind the DNA, they relaxed the region they
# affect...

class TestSuperhelical(TestCase):
    # TODO: Test bindings but superhelical level keeps steady. Here, enzymes may bind and/or unbind but they will
    #  not have effects on the DNA. We need to test for circular and linear. And for effective_size <= size in
    #  some cases.
    #  Test 1: One enzyme binds. We check that global superhelical level doesn't change
    #  Test 2: One enzymbe binds and then it unbinds. We check that the superhelical level do not change.
    #  Test 3: One enzyme binds bare DNA sites. Multiple bindings happen but they do not unbind.
    #           The global superhelical level needs to keep steady.
    #  Test 4: One enzyme that binds bare DNA binds multiple sites, and unbinds. Multiple bindings and unbindings
    #   event happen but the global superhelical level never changes.

    # Tests that updating globals do not change global superhelical density
    def test_update(self):
        my_circuit = Circuit(circuit_circular_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt)

        # This step is very important for defining a new superhelical level
        #for enzyme in my_circuit.enzyme_list:
        #    enzyme.superhelical = superhelical0
        #    enzyme.position = enzyme.position
        #my_circuit.update_twist()
        #my_circuit.update_supercoiling()
        #my_circuit.update_global_twist()
        #my_circuit.update_global_superhelical()

        s0 = my_circuit.superhelical

        for frame in range(1, frames + 1):
            my_circuit.update_twist()
            my_circuit.update_supercoiling()
            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()

            dsigma = abs(s0 - my_circuit.superhelical)
            #print(dsigma)
            self.assertLessEqual(dsigma, error, "Superhelical changed")

    # Test 1: One enzyme binds and it stays bound. We check that global superhelical level doesn't change
    def test_g1(self):
        my_circuit = Circuit(circuit_circular_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt)

        s0 = my_circuit.superhelical

        my_circuit.environmental_list[0].concentration = 10.0

        for frame in range(1, frames + 1):
            # BINDING
            new_enzyme_list = mw.binding_workflow(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                                  my_circuit.rng)
            my_circuit.add_new_enzymes(new_enzyme_list)

            # There are no effects and unbinding
            # UPDATE GLOBALS
            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()

            dsigma = abs(s0 - my_circuit.superhelical)

            if len(new_enzyme_list) > 0:
                print(frame, len(new_enzyme_list), len(my_circuit.enzyme_list), s0, my_circuit.superhelical)
            self.assertLessEqual(dsigma, error, "Superhelical changed")

    #        self.assertEqual(True, False)  # add assertion here

    # TODO: AQUIMEQUEDE. Por aqui esta el error. A veces, cuando unas enzymas bindean el DNA, se dobla el twist,
    #  Y esto cause que se dispare el supercoiling. Piensale y buscale que puede estar pasando.
    #  Test 4: One enzyme that binds bare DNA binds multiple sites, and unbinds. Multiple bindings and unbindings
    #   event happen but the global superhelical level never changes. CIRCULAR CASE
    def test_g4_circular(self):
        my_circuit = Circuit(circuit_circular_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt)

        s0 = my_circuit.superhelical
        twist0 = my_circuit.twist

        my_circuit.environmental_list[2].concentration = 10.0
        # print('START enzyme', len(my_circuit.enzyme_list), 'sigma', my_circuit.superhelical)

        #for frame in range(1, 10):
        for frame in range(1, frames * 3 + 1):
            # BINDING
            new_enzyme_list = mw.binding_workflow(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                                  my_circuit.rng)
            # print(frame, 'BEFORE 1 - enzyme', len(my_circuit.enzyme_list), 'sigma', my_circuit.superhelical, 'twist', my_circuit.twist, sum(enzyme.twist for enzyme in my_circuit.enzyme_list) - my_circuit.enzyme_list[0].twist)
            my_circuit.add_new_enzymes(new_enzyme_list)

            #positions_circuit = [obj.position for obj in my_circuit.enzyme_list]
            #twists_circuit = [obj.twist for obj in my_circuit.enzyme_list]
            #positions_enlist1 = [obj.position for obj in new_enzyme_list]
            #twists_en1 = [obj.twist for obj in new_enzyme_list]

            # There are no effects
            #print(frame, 'BEFORE 2 - enzyme', len(my_circuit.enzyme_list), 'sigma', my_circuit.superhelical, 'twist', my_circuit.twist)

            # UPDATE GLOBALS
            my_circuit.update_global_twist()
            #print(frame, 'BEFORE 3 - enzyme', len(my_circuit.enzyme_list), 'sigma', my_circuit.superhelical, 'twist', my_circuit.twist)
            my_circuit.update_global_superhelical()

            #print(frame, 'BEFORE 4 - enzyme', len(my_circuit.enzyme_list), 'sigma', my_circuit.superhelical, 'twist', my_circuit.twist)

            #if abs(s0 - my_circuit.superhelical) > error:
            #    a=2
            #    s=a+3

            #pcircuit = copy.deepcopy(my_circuit)
            #pnew_enzymes = copy.deepcopy(new_enzyme_list)
            #positions_enlist2 = [obj.position for obj in pnew_enzymes]
            #twists_en2 = [obj.twist for obj in new_enzyme_list]
            #twists_pcircuit = [obj.twist for obj in pcircuit.enzyme_list]
            #positions_pcircuit = [obj.position for obj in pcircuit.enzyme_list]

            if frame == frames:
                print('number of enzymes halfway', len(my_circuit.enzyme_list), 'sigma', my_circuit.superhelical)
                my_circuit.environmental_list[2].concentration = 0.0  # This will turn off the binding

            if frame > frames:
                # UNBINDING
                drop_list_index, drop_list_enzyme = mw.unbinding_workflow(my_circuit.enzyme_list, my_circuit.dt,
                                                                          my_circuit.rng)
                my_circuit.drop_enzymes(drop_list_index)
                my_circuit.add_to_environment(drop_list_enzyme)

                # UPDATE GLOBALS
                my_circuit.update_global_twist()
                my_circuit.update_global_superhelical()

                if len(my_circuit.enzyme_list) == 2:
                    continue

        print(len(my_circuit.enzyme_list), my_circuit.superhelical)
        dsigma = abs(s0 - my_circuit.superhelical)
        print(dsigma)
        self.assertLessEqual(dsigma, error, "Superhelical changed")

    #  Test 4: One enzyme that binds bare DNA binds multiple sites, and unbinds. Multiple bindings and unbindings
    #   event happen but the global superhelical level never changes. LINEAR CASE
    def test_g4_linear(self):
        my_circuit = Circuit(circuit_linear_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt)

        s0 = my_circuit.superhelical

        my_circuit.environmental_list[2].concentration = 10.0

        for frame in range(1, frames * 10 + 1):
            # BINDING
            new_enzyme_list = mw.binding_workflow(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                                  my_circuit.rng)
            my_circuit.add_new_enzymes(new_enzyme_list)

            # There are no effects

            # UPDATE GLOBALS
            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()

            if frame == frames:
                print('number of enzymes halfway', len(my_circuit.enzyme_list), 'sigma', my_circuit.superhelical)
                my_circuit.environmental_list[2].concentration = 0.0  # This will turn off the binding

            if frame > frames:
                # UNBINDING
                drop_list_index, drop_list_enzyme = mw.unbinding_workflow(my_circuit.enzyme_list, my_circuit.dt,
                                                                          my_circuit.rng)
                my_circuit.drop_enzymes(drop_list_index)
                my_circuit.add_to_environment(drop_list_enzyme)

                # UPDATE GLOBALS
                my_circuit.update_global_twist()
                my_circuit.update_global_superhelical()

                if len(my_circuit.enzyme_list) == 2:
                    continue

        print(len(my_circuit.enzyme_list), my_circuit.superhelical)
        dsigma = abs(s0 - my_circuit.superhelical)
        self.assertLessEqual(dsigma, error, "Superhelical changed")
