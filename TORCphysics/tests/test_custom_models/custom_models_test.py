from unittest import TestCase
from TORCphysics import Circuit, Site, Environment, utils
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm
import pandas as pd

# Inputs
circuit_filename = '../test_circuit/circuit.csv'
sites_filename = '../test_circuit/sites.csv'
enzymes_filename = '../test_circuit/enzymes.csv'
environment_filename = '../test_circuit/environment.csv'
output_prefix = ''
frames = 300
series = True
continuation = False
dt = 1
# Define custom models

# Let's define our new enzyme's binding model: CustomEnzymeBinding
class CustomEnzymeBinding(bm.BindingModel):

    # oparams is a dictionary with parameterisation that you can pass to the model
    def __init__(self, filename=None, **oparams):

        super().__init__(filename, **oparams)  # Call the base class constructor

        if not oparams:
            if filename is None:  # If no params are being given
                self.k_on = 0.03  # Rate at which enzyme binds
                self.binding_threshold = -0.01  # Threshold at which the enzyme will bind
            else: # If a file with parameterisation is being given.
                mydata = pd.read_csv(filename)
                if 'k_on' in mydata.columns:
                    self.k_on = mydata['k_on'][0]
                else:
                    raise ValueError('Error, k_on parameter missing in csv file for CustomEnzymeBinding')
                if 'binding_threshold' in mydata.columns:
                    self.binding_threshold = mydata['binding_threshold'][0]
                else:
                    raise ValueError('Error, binding_threshold parameter missing in csv file for CustomEnzymeBinding')
        else:  # If parametersation is being given.
            self.k_on = oparams['k_on']
            self.binding_threshold = oparams['binding_threshold']

        self.oparams = {'k_on': self.k_on, 'binding_threshold': self.binding_threshold}  # Just in case

    # This function calculates the binding probability
    def binding_probability(self, environmental, superhelical, dt):
        # utils is a module that contains useful functions!
        # utils.Poisson_process returns the probability according rate = k_on * concentration
        # and time step dt.
        if superhelical <= self.binding_threshold:
            return utils.Poisson_process(self.k_on * environmental.concentration, dt)
        else:
            return 0.0

    # Lastly, we need this function that tells the model how the rate is modulated as function of supercoiling.
    def rate_modulation(self, superhelical) -> float:
        if superhelical <= self.binding_threshold:
            return self.k_on
        else:
            return 0.0

# We need to initialize the class, then define the calculate_effect() function!
# Let's define our new enzyme's effect model
class CustomEnzymeEffect(em.EffectModel):

    def __init__(self, filename=None, **oparams):

        super().__init__(filename, **oparams)  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.k_cat = 0.001  # Rate at which supercoils are added (supercoils/second)
            else:
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for CustomEnzymeEffect')
        else:
            self.k_cat = oparams['k_cat']

        self.oparams = {'k_cat': self.k_cat}  # Just in case

    # This function calculates the local effect.
    # The model workflow always passes 4 inputs: index = current enzyme index in z_list, current enzyme z, list of enzymes z_list,
    # and time-step dt.
    # Effects are quantified through the effect object em.Effect()
    # Effects are accounted through the change in position (position), twist on the left (twist_left) and twist on the right (twist_right)
    # respect our custom enzyme z.
    def calculate_effect(self, index, z, z_list, dt):
        position = 0.0   # Since the custom enzyme does not move, position=0.0
        twist_left = self.k_cat * dt
        twist_right = self.k_cat * dt
        return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

# We need to initialize the class, then define the unbinding_probability() function!
# Let's define our new enzyme's unbinding model
class CustomEnzymeUnBinding(ubm.UnBindingModel):

    def __init__(self, filename=None, **oparams):

        super().__init__(filename, **oparams)  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.unbinding_threshold = 0.0  # Superhelical level at which the enzyme would unbind.
            else:
                mydata = pd.read_csv(filename)
                if 'unbinding_threshold' in mydata.columns:
                    self.unbinding_threshold = mydata['unbinding_threshold'][0]
                else:
                    raise ValueError('Error, unbinding_threshold parameter missing in csv file for CustomEnzymeUnBinding')
        else:
            self.unbinding_threshold = oparams['unbinding_threshold']

        self.oparams = {'unbinding_threshold': self.unbinding_threshold}  # Just in case

    # This function calculates the unbinding probability.
    # Its inputs are the current enzyme and the time step dt.
    # It returns the unbinding probability; 1 means unbinding and 0.0 will not unbind
    def unbinding_probability(self, enzyme, dt):
        if enzyme.superhelical >= self.unbinding_threshold:
            return 1.0
        else:
            return 0.0

# The actual test

class TestCustomModels(TestCase):

    # We need more tests and better names!
    def test_custom_1(self):

        # Initialize circuit with the initial conditions
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt)

        # We now need to define a binding site for our new enzyme:
        custom_site = Site(site_type='custom_site', name='custom_site_1', start=2500, end=0, k_on=0.03,
                           binding_model=CustomEnzymeBinding())
        #my_circuit.site_list.append(custom_site)

        # Define new environment custom enzyme which can bind its corresponding binding site.
        custom_environment = Environment(e_type='custom_enzyme', name='custom_enzyme_1', site_list=my_circuit.site_list,
                                    concentration=1.0, size=20,
                                    effective_size=10, site_type='custom_site',
                                    effect_model=CustomEnzymeEffect(),
                                    unbinding_model=CustomEnzymeUnBinding())

        my_circuit.add_custom_Site(custom_site)
        my_circuit.add_custom_Environment(custom_environment)
        #my_circuit.environmental_list.append(custom_enzyme)

        # Update lists and twist/supercoiling in the circuit components
        #my_circuit.sort_lists()
        #my_circuit.update_global_twist()
        #my_circuit.update_global_superhelical()
        my_circuit.run()




