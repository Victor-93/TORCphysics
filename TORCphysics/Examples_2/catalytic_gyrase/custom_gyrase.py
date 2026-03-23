import sys
#sys.path.append("/mnt/parscratch/users/username")  # For stanage - replace with you username
from TORCphysics import utils, params
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm
import pandas as pd
import numpy as np
import random

# TODO: The effects need to connect with the environment so it depends on ATP, and that it depends on the force F.
class GyraseCyclesEffect(em.EffectModel):

    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        self.state = 'UNWRAPPED' # The initial state

        # This is just to load some parameters, but won't need it right now.
        if not oparams:
            if filename is None:
                self.k_cat = 10.0
                self.k_wrap = 0.5
                self.k_unwrap = 0.5
                self.k_go = 0.5
                self.k_dwell = 0.5
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
            self.k_dwell = float(oparams['k_dwell'])

            # For making it linear
            #self.a = float(oparams['a'])
            #self.b = float(oparams['b'])

        self.oparams = {'k_cat': self.k_cat, 'k_wrap': self.k_wrap, 'k_unwrap': self.k_unwrap,
                        'k_go': self.k_go, 'k_dwell': self.k_dwell}  # Just in case
                        #'a': self.a, 'b': self.b}

    def calculate_effect(self, index, z, z_list, dt) -> em.Effect:

        position = 0 # This is how much it moves (it doesn't move!)

        # Initial UNWRAPPED state
        # -----------------------------------------------------------------
        # From unwrapped, it can transition to wrap or unbind, but the unbind is decided in the unbind model
        if self.state == 'UNWRAPPED':

            z.name = 'gyrase_' + self.state

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

            z.name = 'gyrase_' + self.state

            # Calculate probabilities as Poisson processes
            p_unwrapped = utils.Poisson_process(self.k_unwrap, dt) # probability of unwrapping
            p_active = utils.Poisson_process(self.k_go, dt)  # probability of transitioning to the active state

            if p_unwrapped + p_active >= 1.0:  # Check that probabilities make sense
                raise ValueError('Error. p_unwrapped + p_active should be less than 1 in GyraseCyclesEffect.')

            # Generate a random number between 0 and 1 to help us decide if it'll transition or not
            random_number = random.random()
            if random_number < p_unwrapped:
                self.state = 'UNWRAPPED'  # TRANSITIONS TO UNWRAPPED
            elif p_unwrapped <= random_number < p_unwrapped + p_active:
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

            z.name = 'gyrase_' + self.state

            # Get superhelical density in the region
            # superhelical = utils.get_superhelical_in_region(z, z_list)

            # Calculate probabilities as Poisson processes
            p_remains = utils.Poisson_process(self.k_dwell, dt)  # probability of unwrapping

            # Generate a random number between 0 and 1 to help us decide if it'll transition or not
            random_number = random.random()
            if random_number < p_remains:
                twist_left = -0.5 * self.k_cat * params.w0 * dt
                twist_right = -0.5 * self.k_cat * params.w0 * dt

                # For the linear model!
                #twist_left = -0.5(self.a + self.b *superhelical* self.k_cat )* params.w0 * dt
                #twist_right = -0.5(self.a + self.b *superhelical* self.k_cat )* params.w0 * dt
            else:
                self.state = 'UNWRAPPED'  # TRANSITIONS TO UNWRAPPED
                twist_left = 0.0
                twist_right = 0.0

            return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

        #return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

# This model is force dependant (stretching force).
# It implements the following kinetic model:
# When Gyrase bounds, it starts  at the UNWRAPPED state.
# In the UNWRAPPED state, it can unbind with k_off, or it can transition to the WRAPPED state,
# which is modelled with a k_wrap rate that is force dependant F and x_wrap.
# x_wrap is the distance to the transition state
# In the WRAPPED state, it can transition to UNWRAPPED according k_unwrap,
# or it can perform the catalytic cycle according a k_go.
# If it performs the catalytic cycle, it returns to the UNWRAPPED state
class GyraseCyclesForce(em.EffectModel):

    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        self.state = 'UNWRAPPED' # The initial state

        # This is just to load some parameters, but won't need it right now.
        if not oparams:
            if filename is None:
                self.k_cat = 20.0
                self.k_wrap = 0.5
                self.k_unwrap = 0.5
                self.k_go = 0.5
                self.force = params.f_stretching
                self.x_wrap = params.kBT_pN_nm * 0.1 # Just guessing
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
            self.force = float(oparams['force'])
            self.x_wrap = float(oparams['x_wrap'])

            # For making it linear
            #self.a = float(oparams['a'])
            #self.b = float(oparams['b'])

        self.oparams = {'k_cat': self.k_cat, 'k_wrap': self.k_wrap, 'k_unwrap': self.k_unwrap,
                        'k_go': self.k_go, 'force': self.force, 'x_wrap':self.x_wrap}  # Just in case
                        #'a': self.a, 'b': self.b}

    def calculate_effect(self, index, z, z_list, dt) -> em.Effect:

        position = 0 # This is how much it moves (it doesn't move!)

        # Initial UNWRAPPED state
        # -----------------------------------------------------------------
        # From unwrapped, it can transition to wrap or unbind, but the unbind is decided in the unbind model
        if self.state == 'UNWRAPPED':

            z.name = 'gyrase_' + self.state
            B = self.force * self.x_wrap / params.kBT_pN_nm
            rate = self.k_wrap * (1 -B)
            probability = utils.P_binding_Nonh_Poisson(rate=rate, dt=dt)

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

            z.name = 'gyrase_' + self.state

            # Calculate probabilities as Poisson processes
            p_unwrapped = utils.Poisson_process(self.k_unwrap, dt) # probability of unwrapping
            p_go = utils.Poisson_process(self.k_go, dt)  # probability of transitioning to the active state

            if p_unwrapped + p_go >= 1.0:  # Check that probabilities make sense
                raise ValueError('Error. p_unwrapped + p_go should be less than 1 in GyraseCyclesForce.')

            # Generate a random number between 0 and 1 to help us decide if it'll transition or not
            random_number = random.random()
            if random_number < p_unwrapped:
                self.state = 'UNWRAPPED'  # TRANSITIONS TO UNWRAPPED
                # No twist is added
                twist_left = 0.0
                twist_right = 0.0
                return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

            elif p_unwrapped <= random_number < p_unwrapped + p_go:
                self.state = 'UNWRAPPED'  # Also transitions to UNWRAPPED
                # But this time, some twist is introduced
                # TODO: This assumes that gyrase is a barrier, which is not true.. We can fix it later.
                twist_left = -0.5 * self.k_cat * params.w0 * dt
                twist_right = -0.5 * self.k_cat * params.w0 * dt
                return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)
            else:
                # Else, do nothing... it stays as WRAPPED

                # The adding twist on the left and the right is 0.0 - this means that it is a barrier
                twist_left = 0.0
                twist_right = 0.0

                return em.Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)



# TODO: Needs additional functions
class GyraseCyclesUnbinding(ubm.UnBindingModel):

    def __init__(self, filename=None, **oparams):

        if not oparams:
            if filename is None:
                self.k_off = 0.1
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_off' in mydata.columns:
                    self.velocity = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_off parameter missing in csv file for GyraseCyclesUnbinding')  # ', filename)
        else:
            self.k_off = float(oparams['k_off'])
        self.oparams = {'k_off': self.k_off}  # Just in case

    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, enzyme, dt) -> float:

        if enzyme.effect_model.state == 'UNWRAPPED':
            rate = self.k_off
            probability = utils.Poisson_process(rate=rate, dt=dt)
        else:
            probability = 0.0
        return probability
