from TORCphysics import utils, params
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm
import pandas as pd
import numpy as np
import random
import sys

# TODO: The effects need to connect with the environment so it depends on ATP, and that it depends on the force F.
class GyraseCyclesEffect(em.EffectModel):

    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        self.state = 'UNWRAPPED' # The initial state

        if not oparams:
            if filename is None:
                self.k_cat = 0.01
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.velocity = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for GyraseCycles')  # ', filename)
        else:
            self.k_cat = float(oparams['k_cat'])
            self.k_wrap = float(oparams['k_wrap'])
            self.k_unwrap = float(oparams['k_unwrap'])
            self.k_go = float(oparams['k_go'])
            self.k_duration = float(oparams['k_duration'])

        self.oparams = {'k_cat': self.k_cat, 'k_wrap': self.k_wrap, 'k_unwrap': self.k_unwrap,
                        'k_go': self.k_go, 'k_duration': self.k_duration}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> em.Effect:

        position = 0 # This is how much it moves (it doesn't move!)

        # Initial UNWRAPPED state
        # -----------------------------------------------------------------
        # From unwrapped, it can transition to wrap or unbind, but the unbind is decided in the unbind model
        if self.state == 'UNWRAPPED':

            rate = self.k_wrap
            probability =  utils.Poisson_process(rate=rate, dt=dt) # Probability modeled as poisson process

            random_number = random.random()
            if random_number <= probability:
                self.state = 'WRAPPED'  # TRANSITIONS TO WRAPPED!

            # Finally, calculate change in twist given the fact that it doesn't form a barrier
            twist_left, twist_right = utils.instant_twist_transfer(z, z_list)

            return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

        # WRAPPED
        # -------------------------------------------------
        # Posibilities, it either transitions to ACTIVE, stays in the WRAPPED state or UNWRAPS.
        if self.state == 'WRAPPED':

            # Calculate probabilities as Poisson processes
            p_unwrapped = utils.Poisson_process(self.k_unwrap, dt) # probability of unwrapping
            p_active = utils.Poisson_process(self.k_go, dt)  # probability of transitioning to the active state

            if p_unwrapped + p_active >= 1.0:  # Check that probabilities make sense
                raise ValueError('Error. p_unwrapped + p_active should be less than 1 in GyraseCyclesEffect.')

            # Generate a random number between 0 and 1 to help us decide if it'll transition or not
            random_number = random.random()
            if random_number < p_unwrapped:
                self.state = 'UNWRAPPED'  # TRANSITIONS TO UNWRAPPED
            elif p_active <= random_number < p_unwrapped + p_active:
                self.state = 'ACTIVE'
            # Else, do nothing... it stays as WRAPPED

            # The adding twist on the left and the right is 0.0 - this means that it is a barrier
            twist_left = 0.0
            twist_right = 0.0

            return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

        # ACTIVE
        # -------------------------------------------------
        # Posibilities: It either keeps acting or stops and returns to the unwrapped state
        if self.state == 'ACTIVE':

            # Calculate probabilities as Poisson processes
            p_remains = utils.Poisson_process(self.k_duration, dt)  # probability of unwrapping

            # Generate a random number between 0 and 1 to help us decide if it'll transition or not
            random_number = random.random()
            if random_number >= p_remains:
                self.state = 'UNWRAPPED'  # TRANSITIONS TO UNWRAPPED

            twist_left = -0.5 * self.k_cat * params.w0 * dt
            twist_right = -0.5 * self.k_cat * params.w0 * dt

            return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


# TODO: Needs additional functions
class GyraseCyclesUnbinding(ubm.UnBindingModel):

    def __init__(self, filename=None, **oparams):
        self.k_off = 0.001  # 0.01  # Rate at which bridge separates
        self.oparams = {'k_off': self.k_off}  # Just in case

    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, enzyme, dt) -> float:

        if enzyme.effect_model.state == 'UNWRAPPED':
            rate = self.k_off
            probability = utils.Poisson_process(rate=rate, dt=dt)
        else:
            probability = 0.0
        return probability
