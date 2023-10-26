from TORCphysics import Circuit, utils
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm
import numpy as np
import random
import sys
from TORCphysics import visualization as vs


# LAC random bridging.
class LacIPoissonBridging(em.EffectModel):

    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        self.k_on = 0.02  # Rate at which bridge forms
        self.k_off = 0.001  # Rate at which bridge separates
        self.state = 'OFF'  # State of the bridge. OFF = No Bridge, ON = Bridge
        self.bridge = False  # True if bridge is formed
        self.bound_with = None  # This is the enzyme that the bridge is bound with
        self.k_cat = 0.01  # Percentage of twist injected (when not bridged)

        # self.oparams = {'velocity': self.velocity, 'gamma': self.gamma}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> em.Effect:

        # TODO: We have to find a way to avoid doing the same bridge calculation for the two enzymes that
        #  form the bridge

        # Everything 0 for now.
        position = 0.0

        rng = np.random.default_rng(random.randrange(sys.maxsize))

        # Bridge is formed
        # -----------------------------------------------------------------
        if self.bridge:

            # Let's calculate the probability of breaking the bridge
            # -----------------------------------------------------------------
            undo_bridge = self.bridge_break(z=z, z_list=z_list, dt=dt, rng=rng)

            # And let's see what will happen to the bridge
            # -----------------------------------------------------------------
            # If it the bridge will undo, then update the state and twist leaks.
            if undo_bridge:
                twist_left, twist_right = self.leak_twist(z=z, z_list=z_list, dt=dt)
                return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

            # If the bridge remains, then nothing happens and not twist leaks.
            else:
                twist_left = 0.0
                twist_right = 0.0
                return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

        # If bridge is not formed
        # -----------------------------------------------------------------
        else:

            # Let's calculate the probability of forming the bridge
            # -----------------------------------------------------------------
            do_bridge = self.bridge_formation(z, z_list, dt, rng)  # This function also updates the state!

            # And let's see what will happen to the bridge
            # -----------------------------------------------------------------
            # If it the bridge will be formed, no twist leaks. The state was already updated previously
            if do_bridge:
                twist_left = 0.0
                twist_right = 0.0
                return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

            # If the bridge does not form, then the state is not updated and twist leaks.
            else:
                twist_left, twist_right = self.leak_twist(z=z, z_list=z_list, dt=dt)
                return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

    def bridge_formation(self, z, z_list, dt, rng):

        # Here, we assume that the bridge is not formed, and we calculate the probability of forming it.

        # First, let's see if there are enzymes of the same type on the enzyme list z_list.
        # -----------------------------------------------------------------
        other_e = [enzyme for enzyme in z_list if enzyme.enzyme_type == z.enzyme_type and enzyme.position != z.position]

        # If there aren't other lacI bound, then don't do anything.
        if len(other_e) == 0:
            do_bridge = False
            return do_bridge

        # If there are other lacI's or candidates for bridging, then let's calculate the probability of
        # forming the bridge
        # -----------------------------------------------------------------
        probability = utils.Poisson_process(rate=self.k_on, dt=dt)

        urandom = rng.uniform()  # we need a random number

        if urandom <= probability:  # and decide if the bridge will form
            do_bridge = True
        else:
            do_bridge = False

        # If the bridge will be formed, select one of the other candidate enzymes other_e to form the bridge
        # and update the state
        if do_bridge:
            # Pick a candidate enzyme with uniform distribution
            random_enzyme = random.choice(other_e)

            # Link and update
            self.bound_with = random_enzyme
            self.bridge = True
            self.state = 'ON'
            # Let's cheat a bit so we can plot the animation and be able to see the site
            # TODO: I think this function is still a bit wrong as sometimes it gives me a double bridge...
            # z.name = z.name+'_bridge'
            z.name = 'lacI_bridge'
            # Then for the other enzyme
            random_enzyme.effect_model.bound_width = z
            random_enzyme.effect_model.bridge = True
            random_enzyme.effect_model.state = 'ON'
            #            random_enzyme.name = random_enzyme.name+'_bridge'
            random_enzyme.name = 'lacI_bridge'

        return do_bridge

    def bridge_break(self, z, z_list, dt, rng):

        # Here, we assume that the bridge is formed, and we calculate the probability of breaking it.

        # Let's calculate the probability of breaking the bridge
        # -----------------------------------------------------------------
        probability = utils.Poisson_process(rate=self.k_off, dt=dt)

        urandom = rng.uniform()  # we need a random number

        if urandom <= probability:  # and decide if the bridge will disappear
            undo_bridge = True
        else:
            undo_bridge = False

        # If the bridge breaks, then update the state
        # -----------------------------------------------------------------
        if undo_bridge:
            # First for the linked enzyme; lets unliked it to this one
            self.bound_with.effect_model.bound_width = None
            self.bound_with.effect_model.bridge = False
            self.bound_with.effect_model.bridge = 'OFF'
            self.bound_with.name = 'lacI'
            # Then unlink and update our enzyme
            self.bound_with = None
            self.bridge = False
            self.state = 'OFF'
            z.name = 'lacI'

        return undo_bridge

    def leak_twist(self, z, z_list, dt):

        # get enzyme on the left of z
        # z_left = utils.get_enzyme_before_position(position=z.position, enzyme_list=z_list)
        # For some reason it fails.... and I have to do the -d. I'll have to check this in the future.
        # If I don't do this, it gives me the incorrect enzyme (it actually gives me the same...)
        d = 0.002
        z_left = [enzyme for enzyme in z_list if enzyme.position <= z.position - d][-1]

        # If the superhelical density on the left is higher than the one on the right, then share a bit...
        if abs(z_left.superhelical) > abs(z.superhelical):
            twist_left = -z_left.twist * self.k_cat * dt
            twist_right = -twist_left
        elif abs(z_left.superhelical) < abs(z.superhelical):
            twist_right = -z.twist * self.k_cat * dt
            twist_left = -twist_right
        # Or if something else happens, don't do anything
        else:
            twist_right = 0.0
            twist_left = 0.0

        return twist_left, twist_right


class LacISimpleUnBinding(ubm.UnBindingModel):

    def __init__(self, filename=None, **oparams):
        self.k_off = 0.001  # 0.01  # Rate at which bridge separates
        self.oparams = {'k_off': self.k_off}  # Just in case

    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, enzyme, dt) -> float:

        # If the bridge is formed, then it can't unbind
        if enzyme.effect_model.bridge:
            probability = 0.0
        # If the bridge isn't formed, then it can unbind
        else:
            probability = utils.Poisson_process(self.k_off, dt)
        return probability


# TODO: El problema de que a veces se queda roja una y la otra verde, puede ser porque a veces una se unbindea cuando
#  se deshac el bridge. Entonces debo de incluir otro binding model, en el que no se pueda unbindear si el bridge
#  esta formado. Lo llamare lacIUnBinding(), que es PoissonProcess pero solo si el bridge no esta formado.


# Let's initialize circuit
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = ''
frames = 1000
series = True
continuation = False
dt = 1.
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
my_circuit.environmental_list[3].effect_model = LacIPoissonBridging()
my_circuit.environmental_list[3].unbinding_model = LacISimpleUnBinding()
# my_circuit.enzyme_list[1].effect_model = LacIPoissonBridging()
my_circuit.run()
my_circuit.print_general_information()

colors_dict = {'tetA': '#d8b41a', 'CDS': 'silver', 'mKalama1': '#0051ff', 'Raspberry': '#e30000'}
output = 'animation'
out_format = '.gif'

vs.create_animation_linear(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, output, out_format,
                           site_type='gene', site_colours=colors_dict)
