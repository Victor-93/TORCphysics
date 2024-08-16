import numpy as np
from TORCphysics import params, utils
import pandas as pd
from abc import ABC, abstractmethod
import sys
import random

# TODO: Decide which of these parameters you need
# All parameters are already in the params module, but I prefer to have them here with more simple names:
v0 = params.v0
w0 = params.w0
gamma = params.gamma


# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains the mathematical functions that compose the effects model.
# Effects can describe the motion of RNAPs as well as twist injection; topoisomerase activity; overall mechanics of
# enzymes bound to DNA, etc...

# ---------------------------------------------------------------------------------------------------------------------
# EFFECT
# ---------------------------------------------------------------------------------------------------------------------
# I thought it'll be easier to describe the effects as an object.
# Because the current effects are taken place at the current enzyme i=index, I use the index to locate the enzyme
# in the enzyme_list which is applying the effect.
# These effects act locally, so they can modify the enzyme's position, and the twist at the neighbouring domains.
class Effect:
    """
    A class used to represent the Effects of bound Enzymes on the DNA molecule.
    These Effects describe local changes on the DNA. These changes include the change in the Enzyme's position,
    and the change in twist on the left/right of the given Enzyme.

    Attributes
    ----------
    index : int
        This is the index that locate the current Enzyme in the list of enzymes 'enzyme_list'.
    position : float
        Parameter that indicates the change in position of the given Enzyme in base-pairs (bp).
        The position of the enzyme after the takes place would be Enzyme.position + position
    twist_left : float
        Parameter that indicates the amount of twist generated on the domain at the left of the given
        Enzyme in radians (rad).
    twist_right : float
        Parameter that indicates the amount of twist generated on the domain at the right of the given
        Enzyme in radians (rad).
    """

    def __init__(self, index, position, twist_left, twist_right):
        """ The constructor of Effect class.

        Parameters
        ----------
        index : int
            This is the index that locate the current Enzyme in the list of enzymes 'enzyme_list'.
        position : float
            Parameter that indicates the change in position of the given Enzyme in base-pairs (bp).
            The position of the enzyme after the takes place would be Enzyme.position + position
        twist_left : float
            Parameter that indicates the amount of twist generated on the domain at the left of the given
            Enzyme in radians (rad).
        twist_right : float
            Parameter that indicates the amount of twist generated on the domain at the right of the given
            Enzyme in radians (rad).
        """

        # I'll save the input filenames just in case
        self.index = index
        self.position = position
        self.twist_left = twist_left
        self.twist_right = twist_right


# ---------------------------------------------------------------------------------------------------------------------
# EFFECT MODELS
# ---------------------------------------------------------------------------------------------------------------------
class EffectModel(ABC):
    """
     The EffectModel abstract class used for defining effect models (subclasses).
     If you need a new model, define it below.
     See how some of the models are defined from this class, so you can make your own and implement it.

     Attributes
     ----------
     filename : str, optional
         Path to the site csv file that parametrises the effect model.
     oparams : dict, optional
         A dictionary containing the parameters used for the effect model.
    """

    def __init__(self, filename=None, continuum=False, **oparams):
        """ The constructor of EffectModel.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the effect model.
        continuum : bool, optional
            Indicates if the actions of the effect model are continuous. Only environmentals should have these type
            of effects, where they don't bind the DNA, but they affect every local region in the molecule continuously.
            For example, you might want the effect of topoisomerases to be continuous.
            Note that if you give an Enzyme that binds the DNA a continuous model, it will only affect it's local
            domain, but the environmental will affect all the local domains.
        oparams : dict, optional
            A dictionary containing the parameters used for the effect model.
        """
        self.filename = filename
        self.continuum = continuum
        self.oparams = oparams

    @abstractmethod
    def calculate_effect(self) -> Effect:
        """ Abstract method for calculating the effect of the Enzyme/molecule.
        This is an essential function for EffectModels as they must be able to calculate the "effect" a given
        Enzyme has on the DNA.

        Returns
        ----------
        effect : Effect
            These functions return an Effect object. This object indicates the enzyme's change in position, and how
            it twists/untwists the DNA on each side for a given timestep.
            Other functions/protocols should then interpret and implement this result.
        """
        pass


# ----------------------
# YOU can define your own models here!
# ----------------------

class RNAPUniform(EffectModel):
    """
     An EffectModel subclass that calculates represents the uniform motion of an RNA Polymerase, while
     injecting positive and negative supercoils (twin domain model).
     This is one of the simplest effect models, where RNAPs can move along the DNA at constant velocity.

     Attributes
     ----------
     velocity : float
        Absolute velocity at which the RNAP moves along the DNA in base-pairs per second (bp/s).
     gamma : float
        Parameter that quantifies the amount of twist generated per base-pair transcribed (rad/bp).
     filename : str, optional
        Path to the site csv file that parametrises the effect model.
     oparams : dict, optional
        A dictionary containing the parameters used for the effect model.
    """

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):
        """ The constructor of the RNAPUniform subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the RNAPUniform effect model; this file should have
            the velocity and gamma parameters
        continuum : bool, optional
            Indicates if the actions of the effect model are continuous. For this model, it is not continuous.
        oparams : dict, optional
            A dictionary containing the parameters used for the effect model. In this case it would be velocity and
            gamma.
        """

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.velocity = params.v0
                self.gamma = params.gamma
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'velocity' in mydata.columns:
                    #  self.velocity = float(mydata['velocity'])
                    self.velocity = mydata['velocity'][0]
                else:
                    raise ValueError('Error, velocity parameter missing in csv file for RNAPUniform')  # ', filename)
                if 'gamma' in mydata.columns:
                    #  self.gamma = float(mydata['gamma'])
                    self.gamma = mydata['gamma'][0]
                else:
                    raise ValueError('Error, gamma parameter missing in csv file for RNAPUniform')  #: ', filename)
        else:
            self.velocity = float(oparams['velocity'])
            self.gamma = float(oparams['gamma'])

        self.oparams = {'velocity': self.velocity, 'gamma': self.gamma}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:
        """ Method for calculating the Effect that the bound RNAP causes on the DNA.

        Parameters
        ----------
        index : int
            Enzyme's index in the list of enzymes "enzyme_list".
        z : Enzyme
            This is the object of the current Enzyme (RNAP) that is moving along the DNA.
        z_list : list
            This is a list of Enzyme objects.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        effect : Effect
            This function returns an Effect object, which indicates the changes in position and local twist that
            the current RNAP caused on the DNA.
        """
        # Get neighbour enzyme
        if z.direction > 0:
            z_n = [e for e in z_list if e.position > z.position][0]  # after - On the right
        if z.direction < 0:
            z_n = [e for e in z_list if e.position < z.position][-1]  # before - On the left
        if z.direction == 0:
            raise ValueError('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')

        # This is if the object moves: simple uniform motion
        position, twist_left, twist_right = uniform_motion(z, dt)

        # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid
        # clashing
        if z.direction > 0:  # Moving to the right
            if z_n.position - (z.position + position) <= 0:
                position = 0.0
                twist_left = 0.0
                twist_right = 0.0
        else:  # Moves to the left
            if z_n.position - (z.position + position) >= 0:
                position = 0.0
                twist_left = 0.0
                twist_right = 0.0

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

class RNAPStagesStall(EffectModel):

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):
        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        # TODO: Vamos a primero visualizar la funcion que queremos, luego la programos aqui.
        #       Pero primero veamos los parametros.
        self.state = 'Closed_complex' # Or Open_complex, or Elongation
        # Closed complex <-> Open complex -> elongation
        # We start in the closed complex state. It can transition to open complex with rate k_oc.
        # It can reverse from open complex to closed with rate k_cc.
        # Then, if in open complex state, it can initiate elongation with rate k_ini. Once the RNAP is in the
        # elongation state, it won't transition to other state until it unbinds. It may stall though.
        if not oparams:
            if filename is None:
                # Closed complex
                self.k_closed = 0.02  # 0.01~100secs Goes from open complex to closed complex

                # Open complex
                self.k_open = 0.05  # Goes from closed complex to opean complex
                self.width = .003   # Width and threshold give the shape to the melting energy
                self.threshold = -.04

                # Initiation (start of elongation)
                self.k_ini = 0.1

                # Elongation params
                self.velocity = params.v0  # Medium velocity
                self.gamma = params.gamma
                self.stall_torque = params.stall_torque
                self.kappa = params.RNAP_kappa

            else:  # There is a file!
                mydata = pd.read_csv(filename)
                # Closed complex
                if 'k_closed' in mydata.columns:
                    self.k_closed = mydata['k_closed'][0]
                else:
                    raise ValueError('Error, k_closed parameter missing in csv file for RNAPStagesStall')

                # Open complex
                if 'k_open' in mydata.columns:
                    self.k_open = mydata['k_open'][0]
                else:
                    raise ValueError('Error, k_open parameter missing in csv file for RNAPStagesStall')
                if 'width' in mydata.columns:
                    self.width = mydata['width'][0]
                else:
                    raise ValueError('Error, width parameter missing in csv file for RNAPStagesStall')
                if 'threshold' in mydata.columns:
                    self.threshold = mydata['threshold'][0]
                else:
                    raise ValueError('Error, threshold parameter missing in csv file for RNAPStagesStall')

                        # Initiation
                if 'k_ini' in mydata.columns:
                    self.k_ini = mydata['k_ini'][0]
                else:
                    raise ValueError('Error, k_ini parameter missing in csv file for RNAPStagesStall')

                # Elongation
                if 'velocity' in mydata.columns:
                    self.velocity = mydata['velocity'][0]
                else:
                    raise ValueError('Error, velocity parameter missing in csv file for RNAPStagesStall')  # ', filename)
                if 'gamma' in mydata.columns:
                    self.gamma = mydata['gamma'][0]
                else:
                    raise ValueError('Error, gamma parameter missing in csv file for RNAPStagesStall')  #: ', filename)
        else:

            # Closed complex
            self.k_closed = float(oparams['k_closed'])

            # Open complex
            self.k_open = float(oparams['k_open'])
            self.width = float(oparams['width'])
            self.threshold = float(oparams['threshold'])

            # Initiation (start of elongation)
            self.k_ini = float(oparams['k_ini'])

            # Elongation
            self.velocity = float(oparams['velocity'])
            self.gamma = float(oparams['gamma'])
            self.stall_torque = float(oparams['stall_torque'])
            self.kappa = float(oparams['kappa'])

        self.oparams = {'k_closed': self.k_closed,
                        'k_open': self.k_open, 'width': self.width, 'threshold': self.threshold,
                        'k_ini': self.k_ini, 'velocity': self.velocity, 'gamma': self.gamma,
                        'kappa': self.kappa, 'stall_torque': self.stall_torque}  # Just in case


    def calculate_effect(self, index, z, z_list, dt) -> Effect:

        # TODO: Just figure out how you treat the twist leak, does the RNAP always forms a barrier in all stages?
        #  I was thinking that maybe in the closed complex it leaks the twist, but in the open it doesn't.
        # TODO: Create example of the model.
        # TODO: Twist leak: In the closed_complex state, the RNAP acts as a barrier? Or does it leaks superhelical density?
        #  Maybe we should assume that it leaks with certain rate?
        # TODO: Twist leak2: In the open_complex state, it starts acting as a barrier?
        # Everything 0 for now
        position = 0.0
        twist_left = 0.0
        twist_right = 0.0
        # Closed complex
        # -------------------------------------------------
        # If we are in the closed_complex state, then it can either transition to the open complex state
        # or unbind. But the unbinding_model is in charge of the unbinding condition. So here, we just test if
        # it transitions to the open complex state.
        # Superhelical density for openning the DNA. It is the superhelical density acting on the region
        if self.state == 'Closed_complex':

            # Get superhelical density in the region
            superhelical_region, twist_region = utils.get_superhelical_in_region(z, z_list)
            # Calculate the opening function U (it is not actually the energy because it is scaled with an effective
            # energy)
            U = utils.opening_function(superhelical_region, self.threshold, self.width)

            rate = self.k_open * np.exp(-U)
            probability =  utils.P_binding_Nonh_Poisson(rate=rate, dt=dt)

            # TODO: We have a rng problem! I define it outside this funciton, in the main code (circuit).
            #  Maybe I should pass the rng to here as well, so we always use the same seed?
            # Generate a random number between 0 and 1 to help us decide if it'll form the open complex
            random_number = random.random()
            if random_number <= probability:
                self.state = 'Open_complex'  # TRANSITIONS TO OPEN COMPLEX!

            return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

        # Open complex
        # -------------------------------------------------
        # Posibilities, it either transitions to closed complex, stays as an open complex or transitions to elongation.
        if self.state == 'Open_complex':
            # Calculate probabilities as Poisson processes
            p_closed = utils.Poisson_process(self.k_closed, dt) # probability of forming closed complex
            p_init = utils.Poisson_process(self.k_ini, dt)  # probability of elongation initiation
            # p_open = 1 - (p_closed + p_init)  # Probability of staying as open complex

            if p_closed + p_init >= 1.0:  # Check that probabilities make sense
                raise ValueError('Error. p_closed + p_init should be less than 1 in RNAPStagesStall.')

            # Generate a random number between 0 and 1 to help us decide if it'll transition or not
            random_number = random.random()
            if random_number < p_closed:
                self.state = 'Closed_complex'  # TRANSITIONS TO CLOSED COMPLEX
            elif p_closed <= random_number < p_closed + p_init:
                self.state = 'Elongation'
            # Else, do nothing... it stays as an open complex

            return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

        # Elongation
        # -------------------------------------------------
        if self.state == 'Elongation':
            # Elongation
            position, twist_left, twist_right = utils.RNAP_stall_mec_model(z, z_list, dt)
            return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

        if self.state not in ['Closed_complex', 'Open_complex', 'Elongation']:
            raise ValueError('Error. Unknown state in RNAPStagesStall.')


# TODO: Document the RNAPStall model. It is a model with velocity but no torques.
# TODO: Test this function
class RNAPStall(EffectModel):

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):
        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.velocity = params.v0  # Medium velocity
                self.gamma = params.gamma
                self.stall_torque = params.stall_torque
                self.kappa = params.RNAP_kappa
            else:  # There is a file!
                # TODO: Do it for the other new parameters
                mydata = pd.read_csv(filename)
                if 'velocity' in mydata.columns:
                    #  self.velocity = float(mydata['velocity'])
                    self.velocity = mydata['velocity'][0]
                else:
                    raise ValueError('Error, velocity parameter missing in csv file for RNAPUniform')  # ', filename)
                if 'gamma' in mydata.columns:
                    #  self.gamma = float(mydata['gamma'])
                    self.gamma = mydata['gamma'][0]
                else:
                    raise ValueError('Error, gamma parameter missing in csv file for RNAPUniform')  #: ', filename)
        else:
            self.velocity = float(oparams['velocity'])
            self.gamma = float(oparams['gamma'])
            self.stall_torque = float(oparams['stall_torque'])
            self.kappa = float(oparams['kappa'])

        self.oparams = {'velocity': self.velocity, 'gamma': self.gamma,
                        'kappa': self.kappa, 'stall_torque': self.stall_torque}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:

        # if z.direction == 0:
        #    raise ValueError('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')

        # Get enzymes on the right and left
        z_right = [e for e in z_list if e.position > z.position][0]  # after - On the right
        z_left = [e for e in z_list if e.position < z.position][-1]  # before - On the left

        # First, we need to calculate the Torque acting on our RNAP.
        # Calculate torques and determine if the RNAP will stall
        torque_right = Marko_torque(z.superhelical)  # Torque on the right
        torque_left = Marko_torque(z_left.superhelical)  # Torque on the left

        if z.direction > 0:  # Moving to the right
            torque = torque_right - torque_left
        else:  # Moving to the left
            torque = torque_left - torque_right

        velocity = velocity_2022SevierBioJ(z, torque)

        position = z.direction * velocity * dt

        # Injects twist: denatures w = gamma*velocity*dt base-pairs
        twist_left = -z.direction * z.effect_model.gamma * params.w0 * velocity * dt
        twist_right = z.direction * z.effect_model.gamma * params.w0 * velocity * dt

        # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid
        # clashing
        if z.direction > 0:  # Moving to the right
            if z_right.position - (z.position + position) <= 0:
                position = 0.0
                twist_left = 0.0
                twist_right = 0.0
        else:  # Moves to the left
            if z_left.position - (z.position + position) >= 0:
                position = 0.0
                twist_left = 0.0
                twist_right = 0.0

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


# TODO: Check that effect distributes supercoiling on both side (left and right)
class TopoIUniform(EffectModel):
    """
     An EffectModel subclass that represents the uniform effect that topoisomerase I have on the DNA.
     In this model, bound enzymes inject supercoils/twist to the left and right uniformly, that is, that the
     amounts of supercoils injected are independent of the local supercoiling density.
     For each timestep, supercoils will be injected constantly until the enzyme unbinds.

     Attributes
     ----------
     k_cat : float
        Catalysis rate at which supercoils are being removed per second (bp/sec).
     filename : str, optional
        Path to the site csv file that parametrises the effect model.
     oparams : dict, optional
        A dictionary containing the parameters used for the effect model.
    """

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):
        """ The constructor of the TopoIPUniform subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the TopoIPUniform effect model; this file should have
            the k_cat parameter
        continuum : bool, optional
            Indicates if the actions of the effect model are continuous. For this model, it is not continuous.
        oparams : dict, optional
            A dictionary containing the parameters used for the effect model. In this case it would be k_cat.
        """

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.k_cat = params.topoI_uniform_k_cat
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for TopoIUniform')
        else:
            self.k_cat = float(oparams['k_cat'])

        self.oparams = {'k_cat': self.k_cat}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:
        """ Method for calculating the simple and uniform Effect that the bound Topoisomerase I cause on the DNA.

        Parameters
        ----------
        index : int
            Enzyme's index in the list of enzymes "enzyme_list".
        z : Enzyme
            This is the object of the current Enzyme (RNAP) that is moving along the DNA.
        z_list : list
            This is a list of Enzyme objects.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        effect : Effect
            This function returns an Effect object, which indicates the changes in position and local twist that
            the current Topo I caused on the DNA.
        """

        position, twist_left, twist_right = topoisomerase_supercoiling_injection(self.k_cat, dt)

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


# TODO: Check that effect distributes supercoiling on both side (left and right)
class GyraseUniform(EffectModel):
    """
     An EffectModel subclass that represents the uniform effect that gyrase have on the DNA.
     In this model, bound enzymes inject supercoils/twist to the left and right uniformly, that is, that the
     amounts of supercoils injected are independent of the local supercoiling density.
     For each timestep, supercoils will be injected constantly until the enzyme unbinds.

     Attributes
     ----------
     k_cat : float
        Catalysis rate at which supercoils are being removed per second (bp/sec).
     filename : str, optional
        Path to the site csv file that parametrises the effect model.
     oparams : dict, optional
        A dictionary containing the parameters used for the effect model.
    """

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):
        """ The constructor of the GyrasePUniform subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the GyraseUniform effect model; this file should have
            the k_cat parameter.
        continuum : bool, optional
            Indicates if the actions of the effect model are continuous. For this model, it is not continuous.
        oparams : dict, optional
            A dictionary containing the parameters used for the effect model. In this case it would be k_cat.
        """

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.k_cat = params.gyra_uniform_k_cat
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for GyraseUniform')
        else:
            self.k_cat = float(oparams['k_cat'])

        self.oparams = {'k_cat': self.k_cat}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:
        """ Method for calculating the simple and uniform Effect that the bound Gyrase cause on the DNA.

        Parameters
        ----------
        index : int
            Enzyme's index in the list of enzymes "enzyme_list".
        z : Enzyme
            This is the object of the current Enzyme (RNAP) that is moving along the DNA.
        z_list : list
            This is a list of Enzyme objects.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        effect : Effect
            This function returns an Effect object, which indicates the changes in position and local twist that
            the current Gyrase caused on the DNA.
        """

        position, twist_left, twist_right = topoisomerase_supercoiling_injection(self.k_cat, dt)

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


class GyraseEnergyPack(EffectModel):
    """
     An EffectModel subclass that represents the effect that gyrase have on the DNA.
     In this model, bound gyrase inject supercoils/twist to the left and right as a function of the free energy
     associated with the local superhelical density. In this model, gyrase can only introduce a finite amount
     of negative supercoils to the DNA, that is, bound gyrases have energy packages.
     When the energy of these packages are depleted, gyrase will not be able to introduce supercoils.
     For each timestep, supercoils will be injected constantly until the enzyme unbinds.

     Attributes
     ----------
     k_cat : float
        Catalysis rate at which supercoils are being introduced per second (supercoils/sec = 1/sec).
        * Note: superhelical density is dimensionless
     filename : str, optional
        Path to the site csv file that parametrises the effect model.
     oparams : dict, optional
        A dictionary containing the parameters used for the effect model.
    """

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):
        """ The constructor of the GyrasePUniform subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the GyraseUniform effect model; this file should have
            the k_cat parameter.
        continuum : bool, optional
            Indicates if the actions of the effect model are continuous. For this model, it is not continuous.
        oparams : dict, optional
            A dictionary containing the parameters used for the effect model. In this case it would be k_cat.
        """

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.k_cat = 0.0001  # TODO: determine params.gyra_energypack_k_cat
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for GyraseEnergyPack')
        else:
            self.k_cat = float(oparams['k_cat'])

        self.energy_pack =  params.ATP_hydrolysis  # The amount of free energy it can pay

        self.oparams = {'k_cat': self.k_cat}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:
        """ Method for calculating the Effect that the bound Gyrase cause on the DNA.

        Parameters
        ----------
        index : int
            Enzyme's index in the list of enzymes "enzyme_list".
        z : Enzyme
            This is the object of the current Enzyme (RNAP) that is moving along the DNA.
        z_list : list
            This is a list of Enzyme objects.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        effect : Effect
            This function returns an Effect object, which indicates the changes in position and local twist that
            the current Gyrase caused on the DNA.
        """

        position = 0.0  # Gyrase does not move

        # If we have energy, then spend some in injecting supercoils.
        if self.energy_pack > 0.0:
            # TODO: Think if you can get put this process in the utils file
            # Let's obtain the superhelical density within the domain
            z_before = utils.get_enzyme_before_position(position=z.position - 5, enzyme_list=z_list)  # Get the enzyme before
            z_after = utils.get_enzyme_after_position(position=z.position + 5, enzyme_list=z_list)  # Get the enzyme after z
            total_twist = z_before.twist + z.twist  # Total twist in the region.
            # This is the total superhelical density of a region; enzyme X does not block supercoils  O______X_____O
            superhelical = total_twist / (params.w0 * (z_after.position - z_before.position))
            length = utils.calculate_length(z_before, z_after)

            G_sigma = utils.superhelical_free_energy(superhelical, length)

            # ------
            twist_left = -0.5 * self.k_cat * np.exp(-G_sigma) * dt * params.w0 * length
            supercoils_injected = -self.k_cat * np.exp(-G_sigma) * dt

            # ------
            energy_spent = utils.change_supercoiling_free_energy(superhelical,
                                                                 superhelical + supercoils_injected, length)
            self.update_energy_pack(abs(energy_spent))
        # If we don't have energy, then do nothing.
        else:
            twist_left = 0.0

        twist_right = twist_left
        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

    # If we spent energy, then the stored energy (energy_pack) decreases
    def update_energy_pack(self, energy_spent):
        self.energy_pack = self.energy_pack - energy_spent


# TODO: Check that effect distributes supercoiling on both side (left and right)
# TODO: Comment and fix
class TopoisomeraseLinearEffect(EffectModel):

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        # TODO: Check correct parametrization
        if not oparams:
            if filename is None:
                self.k_cat = params.gyra_uniform_k_cat
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                    self.sigma0 = mydata['sigma0'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for GyraseUniform')
        else:
            self.k_cat = float(oparams['k_cat'])
            self.sigma0 = float(oparams['sigma0'])

        self.oparams = {'k_cat': self.k_cat, 'sigma0': self.sigma0}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:

        position = 0.0
        z_b = utils.get_enzyme_before_position(position=z.position - 10, enzyme_list=z_list)  # Get the enzyme before
        z_a = utils.get_enzyme_after_position(position=z.position + 10, enzyme_list=z_list)  # Get the enzyme after z
        total_twist = z_b.twist + z.twist  # Total twist in the region.
        # This is the total superhelical density of a region; enzyme X does not block supercoils  O______X_____O
        superhelical = total_twist / (params.w0 * (z_a.position - z_b.position))
        twist_left = 0.5 * self.k_cat * params.w0 * dt * (self.sigma0 - superhelical)
        twist_right = twist_left

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


# TODO: Make the random parameter an input.
# TODO: Check that effect distributes supercoiling on both side (left and right)
class TopoisomeraseLinearRandEffect(EffectModel):

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        # TODO: Check correct parametrization
        if not oparams:
            if filename is None:
                self.k_cat = params.gyra_uniform_k_cat
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for GyraseUniform')
        else:
            self.k_cat = float(oparams['k_cat'])
            self.sigma0 = float(oparams['sigma0'])

        self.oparams = {'k_cat': self.k_cat, 'sigma0': self.sigma0}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:

        position = 0.0
        z_b = utils.get_enzyme_before_position(position=z.position - 10, enzyme_list=z_list)  # Get the enzyme before
        z_a = utils.get_enzyme_after_position(position=z.position + 10, enzyme_list=z_list)  # Get the enzyme after z
        total_twist = z_b.twist + z.twist  # Total twist in the region.
        superhelical = total_twist / (params.w0 * (z_a.position - z_b.position))

        random_addition = np.random.uniform(-0.1, 0.1)  # This is a random variation of supercoils introduced

        twist_left = 0.5 * self.k_cat * params.w0 * dt * (self.sigma0 - superhelical + random_addition)
        twist_right = twist_left

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


# TODO: Comment and fix
# The function is : delta twist = kcat * w0 * (sigma0 - sigma) * dt
# sigma0 marks the superhelical density in which the torque of DNA is too strong that it restores the superhelicity
# in other words, gyrase cannot keep acting at that superhelicity
class GyraseLinear(EffectModel):

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        # TODO: Check correct parametrization
        if not oparams:
            if filename is None:
                self.k_cat = params.gyra_uniform_k_cat
                self.sigma0 = -0.2  # TODO: I don't know why like this
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                    self.sigma0 = mydata['sigma0'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for GyraseUniform')
        else:
            self.k_cat = float(oparams['k_cat'])
            self.sigma0 = float(oparams['sigma0'])

        self.oparams = {'k_cat': self.k_cat, 'sigma0': self.sigma0}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:

        position = 0.0
        z_b = utils.get_enzyme_before_position(position=z.position - 10, enzyme_list=z_list)  # Get the enzyme before
        z_a = utils.get_enzyme_after_position(position=z.position + 10, enzyme_list=z_list)  # Get the enzyme after z
        total_twist = z_b.twist + z.twist  # Total twist in the region.
        # This is the total superhelical density of a region; enzyme X does not block supercoils  O______X_____O
        superhelical = total_twist / (params.w0 * (z_a.position - z_b.position))

        # Total twist added by the molecule
        twist_added = self.k_cat * params.w0 * dt * (self.sigma0 - superhelical)

        # Total superhelicity in the region after the effect
        total_superhelicity = (total_twist + twist_added) / (params.w0 * (z_a.position - z_b.position))

        # The twist added is distributed so we have the same total superhelicity at both sides after the effect.
        # delta twist = total_superhelicity * w0 * region_length - twist_before_effect
        twist_left = total_superhelicity * params.w0 * abs(z.position - z_b.position) - z_b.twist
        twist_right = total_superhelicity * params.w0 * abs(z_a.position - z.position) - z.twist

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


    def twist_added(self, superhelical, dt):
        # Total twist added by the molecule
        twist_added = self.k_cat * params.w0 * dt * (self.sigma0 - superhelical)
        return twist_added

# TODO: Comment and fix
# The function is : delta twist = -kcat * w0 * sigma * dt
class TopoILinear(EffectModel):

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        # TODO: Check correct parametrization
        if not oparams:
            if filename is None:
                self.k_cat = params.topoI_uniform_k_cat
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for TopoILinear')
        else:
            self.k_cat = float(oparams['k_cat'])

        self.oparams = {'k_cat': self.k_cat}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:

        position = 0.0
        z_b = utils.get_enzyme_before_position(position=z.position - 10, enzyme_list=z_list)  # Get the enzyme before
        z_a = utils.get_enzyme_after_position(position=z.position + 10, enzyme_list=z_list)  # Get the enzyme after z
        total_twist = z_b.twist + z.twist  # Total twist in the region.
        # This is the total superhelical density of a region; enzyme X does not block supercoils  O______X_____O
        superhelical = total_twist / (params.w0 * (z_a.position - z_b.position))

        # Total twist added by the molecule
        twist_added = -self.k_cat * params.w0 * dt * superhelical

        # Total superhelicity in the region after the effect
        total_superhelicity = (total_twist + twist_added) / (params.w0 * (z_a.position - z_b.position))

        # The twist added is distributed so we have the same total superhelicity at both sides after the effect.
        # delta twist = total_superhelicity * w0 * region_length - twist_before_effect
        twist_left = total_superhelicity * params.w0 * abs(z.position - z_b.position) - z_b.twist
        twist_right = total_superhelicity * params.w0 * abs(z_a.position - z.position) - z.twist

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

    def twist_added(self, superhelical, dt):
        # Total twist added by the molecule
        twist_added = -self.k_cat * params.w0 * dt * superhelical
        return twist_added

# TODO: Comment and fix
# The function is : delta twist = angular_vel0 * dt - k_DNA * w0 * sigma * dt*dt
# k_DNA = acceleration at which supercoils are being removed (bp/sec^2) -
#   This parameter is related with the torque of DNA
class TopoISimpleTorque(EffectModel):

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=False, **oparams):

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        # TODO: Check correct parametrization
        if not oparams:
            if filename is None:
                self.k_DNA = 5.0
                self.drag_constant = 0.0
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_DNA' in mydata.columns:
                    self.k_DNA = mydata['k_DNA'][0]
                else:
                    raise ValueError('Error, k_DNA parameter missing in csv file for TopoISimpleTorque')
        else:
            self.k_DNA = float(oparams['k_DNA'])
            self.drag_constant = float(oparams['drag_constant'])

        self.angular_velocity = 0.0  # Angular velocity
        self.oparams = {'k_DNA': self.k_DNA, 'drag_constant': self.drag_constant}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:

        position = 0.0
        z_b = utils.get_enzyme_before_position(position=z.position - 10, enzyme_list=z_list)  # Get the enzyme before
        z_a = utils.get_enzyme_after_position(position=z.position + 10, enzyme_list=z_list)  # Get the enzyme after z
        total_twist = z_b.twist + z.twist  # Total twist in the region.
        # This is the total superhelical density of a region; enzyme X does not block supercoils  O______X_____O
        superhelical = total_twist / (params.w0 * (z_a.position - z_b.position))

        # Total twist added by the molecule
#        twist_added = self.angular_velocity * dt - self.k_DNA * params.w0 * superhelical * dt * dt
        twist_added = (self.angular_velocity * dt - self.k_DNA * params.w0 * superhelical * dt * dt
                       - self.drag_constant * self.angular_velocity * dt)

        # Total superhelicity in the region after the effect
        total_superhelicity = (total_twist + twist_added) / (params.w0 * (z_a.position - z_b.position))

        # The twist added is distributed so we have the same total superhelicity at both sides after the effect.
        # delta twist = total_superhelicity * w0 * region_length - twist_before_effect
        twist_left = total_superhelicity * params.w0 * abs(z.position - z_b.position) - z_b.twist
        twist_right = total_superhelicity * params.w0 * abs(z_a.position - z.position) - z.twist

        self.update_angular_velocity(twist_added, dt)

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)

    # delta_twist = change in twist
    def update_angular_velocity(self, delta_twist, dt):
        self.angular_velocity = delta_twist/dt


class TopoIContinuum(EffectModel):
    """
     An EffectModel subclass that calculates represents the continuum effect of Topoisomerase I, on the DNA.
     This model affects every region on the DNA continuously. These effects are represented by a sigmoid curve.
     This model is compatible with the Houdaigui et al. 2019 model.

     The amount of supercoils removed is calculated by:
     supercoils_removed = concentration * k_cat * dt / (1 + exp( (supercoiling - threshold)/width)


     Attributes
     ----------
     k_cat : float
        Catalysis rate at which supercoils are being removed per second (1/nM*s).
     threshold : float
        The threshold of the sigmoid curve. This parameter is dimensionless.
     width : float
        The width of the sigmoid curve. This parameter is dimensionless.
     filename : str, optional
        Path to the site csv file that parametrises the effect model.
     oparams : dict, optional
        A dictionary containing the parameters used for the effect model.
    """

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=True, **oparams):
        """ The constructor of the RNAPUniform subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the TopoIContinuum effect model; this file should have
            the k_cat, threshold and width parameters.
        continuum : bool, optional
            Indicates if the actions of the effect model are continuous. For this model, it is!
        oparams : dict, optional
            A dictionary containing the parameters used for the effect model. In this case, it would be k_cat,
            threshold and width
        """

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.k_cat = params.topo_sam_kcat
                self.threshold = params.topo_sam_threshold
                self.width = params.topo_sam_width
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for TopoIContinuum')
                if 'threshold' in mydata.columns:
                    self.threshold = mydata['threshold'][0]
                else:
                    raise ValueError('Error, threshold parameter missing in csv file for TopoIContinuum')
                if 'width' in mydata.columns:
                    self.width = mydata['width'][0]
                else:
                    raise ValueError('Error, width parameter missing in csv file for TopoIContinuum')
        else:
            self.k_cat = float(oparams['k_cat'])
            self.threshold = float(oparams['threshold'])
            self.width = float(oparams['width'])

        self.oparams = {'k_cat': self.k_cat, 'threshold': self.threshold, 'width': self.width}  # Just in case

    def calculate_effect(self, concentration, index, z, z_list, dt) -> Effect:

        """ Method for calculating the Effect that continuum action of TopoI causes on the DNA.

         Parameters
         ----------
         concentration : float
             Enzyme concentration in the environment.
         index : int
             Enzyme's index in the list of enzymes "enzyme_list".
         z : Enzyme
             This is the object of the current Enzyme (RNAP) that is moving along the DNA.
         z_list : list
             This is a list of Enzyme objects.
         dt : float
             Timestep in seconds (s).

         Returns
         ----------
         effect : Effect
             This function returns an Effect object, which indicates the changes in position and local twist that
             TopoI caused on the DNA.
         """

        # Calculates the amount of coils removed by topoisomerase I activity.
        # This function only depends on the supercoiling density (sigma)
        # I took this function from Sam Meyer's paper (2019)
        # the function has the form of (concentration*sigmoid)*rate*dt
        z_n = [e for e in z_list if e.position > z.position][0]  # Enzyme on the right
        a = concentration * self.k_cat * dt
        try:
            b = 1 + np.exp((z.superhelical - self.threshold) / self.width)
            supercoiling_removed = a / b
        except OverflowError as oe:
            supercoiling_removed = 0.0

        twist_right = utils.calculate_twist_from_sigma(z, z_n, supercoiling_removed)
        return Effect(index=index, position=0.0, twist_left=0.0, twist_right=twist_right)


# TODO: Check if it is easier to find the next neighbour? z_n? Maybe a function that can speed things up
class GyraseContinuum(EffectModel):
    """
     An EffectModel subclass that calculates represents the continuum effect of Gyrase on the DNA.
     This model affects every region on the DNA continuously. These effects are represented by a sigmoid curve.
     This model is compatible with the Houdaigui et al. 2019 model.

     The amount of supercoils removed is calculated by:
     supercoils_removed = concentration * k_cat * dt / (1 + exp( (supercoiling - threshold)/width)

     Attributes
     ----------
     k_cat : float
        Catalysis rate at which supercoils are being removed per second (1/nM*s).
     threshold : float
        The threshold of the sigmoid curve. This parameter is dimensionless.
     width : float
        The width of the sigmoid curve. This parameter is dimensionless.
     filename : str, optional
        Path to the site csv file that parametrises the effect model.
     oparams : dict, optional
        A dictionary containing the parameters used for the effect model.
    """

    # def __init__(self, name, filename):
    def __init__(self, filename=None, continuum=True, **oparams):
        """ The constructor of the RNAPUniform subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the TopoIContinuum effect model; this file should have
            the k_cat, threshold and width parameters.
        continuum : bool, optional
            Indicates if the actions of the effect model are continuous. For this model, it is!
        oparams : dict, optional
            A dictionary containing the parameters used for the effect model. In this case, it would be k_cat,
            threshold and width
        """

        super().__init__(filename, continuum, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.k_cat = params.gyra_sam_kcat
                self.threshold = params.gyra_sam_threshold
                self.width = params.gyra_sam_width
            else:  # There is a file!
                mydata = pd.read_csv(filename)
                if 'k_cat' in mydata.columns:
                    self.k_cat = mydata['k_cat'][0]
                else:
                    raise ValueError('Error, k_cat parameter missing in csv file for GyraseContinuum')
                if 'threshold' in mydata.columns:
                    self.threshold = mydata['threshold'][0]
                else:
                    raise ValueError('Error, threshold parameter missing in csv file for GyraseContinuum')
                if 'width' in mydata.columns:
                    self.width = mydata['width'][0]
                else:
                    raise ValueError('Error, width parameter missing in csv file for GyraseContinuum')
        else:
            self.k_cat = float(oparams['k_cat'])
            self.threshold = float(oparams['threshold'])
            self.width = float(oparams['width'])

        self.oparams = {'k_cat': self.k_cat, 'threshold': self.threshold, 'width': self.width}  # Just in case

    def calculate_effect(self, concentration, index, z, z_list, dt) -> Effect:

        """ Method for calculating the Effect that continuum action of Gyrase causes on the DNA.

         Parameters
         ----------
         concentration : float
             Enzyme concentration in the environment.
         index : int
             Enzyme's index in the list of enzymes "enzyme_list".
         z : Enzyme
             This is the object of the current Enzyme (RNAP) that is moving along the DNA.
         z_list : list
             This is a list of Enzyme objects.
         dt : float
             Timestep in seconds (s).

         Returns
         ----------
         effect : Effect
             This function returns an Effect object, which indicates the changes in position and local twist that
             Gyrase caused on the DNA.
         """

        # Calculates the amount of coils removed by gyrase activity.
        # This function only depends on the supercoiling density (sigma)
        # I took this function from Sam Meyer's paper (2019)
        # the function has the form of (concentration*sigmoid)*rate*dt

        z_n = [e for e in z_list if e.position > z.position][0]  # Enzyme on the right
        a = concentration * self.k_cat * dt
        try:
            b = 1 + np.exp(-(z.superhelical - self.threshold) / self.width)
            supercoiling_removed = -a / b
        except OverflowError as oe:
            supercoiling_removed = 0.0
        twist_right = utils.calculate_twist_from_sigma(z, z_n, supercoiling_removed)
        return Effect(index=index, position=0.0, twist_left=0.0, twist_right=twist_right)


# ---------------------------------------------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

# According inputs, loads the effect model, name and its params. This function is used in environment and enzyme.
# This function calls assign_effect_model
def get_effect_model(name, e_model, model_name, oparams_file, oparams):
    """ This function loads the EffectModel to implement according the provided inputs.
    This function is used for Environments and Enzymes. So this function is implemented by those two classes.

    Parameters
    ----------
    name : str
        Name of the environmental or enzyme.
    e_model : EffectModel or None
        An EffectModel or None.
    model_name : str
        Name of the model to use, e.g. 'RNAPUniform'
    oparams_file : str, optional
        Path to the csv file containing the parametrisation of the EffectModel to use.
    oparams : dict, optional
        A dictionary containing the parameters used for the effect model. In the case of RNAPUniform, it would be
        velocity and gamma.

    Returns
    ----------
    effect_model : EffectModel or None
        The EffectModel to implement for the Enzyme/Environment. If no EffectModel could be determined, this variable
        will be None.
    effect_model_name: str or None
        Name of the EffectModel to use. It is the same as effect_model.__class__.__name__
        If the EffectModel was not determined, then this variable is None.
    effect_oparams_file: str or None
        Path to the csv file containing the parametrisation of the EffectModel. None if file was not given.
    effect_model_oparams : dict or None
        Dictionary with the parametrisation of the EffectModel. None will be returned if the EffectModel could not
        be determined.
    """

    # If no model is given
    if e_model is None:

        # No model is given, not even a name, so there's NO effect model
        if model_name is None:
            e_model = None
            model_name = None
            oparams_file = None
            oparams = None

        # Model indicated by name
        else:
            # Loads effect model.
            # If oparams dict is given, those will be assigned to the model -> This is priority over oparams_file
            # If oparams_file is given, parameters will be read from file, in case of no oparams dict
            # If no oparams file/dict are given, default values will be used.

            # A dictionary of parameters is given so that's priority
            if isinstance(oparams, dict):
                e_model = assign_effect_model(model_name, **oparams)
            # No dictionary was given
            else:
                # If no oparams_file is given, then DEFAULT values are used.
                if oparams_file is None:
                    e_model = assign_effect_model(model_name)
                # If an oparams_file is given, then those are loaded
                else:
                    e_model = assign_effect_model(model_name, oparams_file=oparams_file)

                oparams = e_model.oparams  # To make them match

    # An actual model was given
    else:

        #  Let's check if it's actually an effect model - The model should already have the oparams
        if isinstance(e_model, EffectModel):
            #  Then, some variables are fixed.
            model_name = e_model.__class__.__name__
            oparams = e_model.oparams
            oparams_file = None

        else:
            print('Warning, effect model given is not a class for environmental/enzyme ', name)
            e_model = None
            model_name = None
            oparams_file = None
            oparams = None

    effect_model = e_model
    effect_model_name = model_name
    effect_oparams_file = oparams_file
    effect_model_oparams = oparams

    return effect_model, effect_model_name, effect_oparams_file, effect_model_oparams


# Add your models into this function so it the code can recognise it
def assign_effect_model(model_name, oparams_file=None, **oparams):
    """ This function decides the EffectModel to use according the provided inputs.

    Parameters
    ----------

    model_name : str
        Name of the EffectModel to use. e,g, RNAPUniform.
    oparams_file : str, optional
        Path to the csv file containing the parametrisation of the EffectModel to use.
    oparams : dict, optional
        A dictionary containing the parameters used for the effect model. In the case of RNAPUniform, it would be
        velocity and gamma.

    Returns
    ----------
    my_model : EffectModel
        A EffectModel object that describes the effect mechanism of the given Enzyme.
    """
    if model_name == 'RNAPUniform':
        my_model = RNAPUniform(filename=oparams_file, **oparams)
    elif model_name == 'RNAPStall':
        my_model = RNAPStall(filename=oparams_file, **oparams)
    elif model_name == 'TopoIUniform':
        my_model = TopoIUniform(filename=oparams_file, **oparams)
    elif model_name == 'GyraseUniform':
        my_model = GyraseUniform(filename=oparams_file, **oparams)
    elif model_name == 'TopoisomeraseLinearEffect':
        my_model = TopoisomeraseLinearEffect(filename=oparams_file, **oparams)
    elif model_name == 'TopoisomeraseLinearRandEffect':
        my_model = TopoisomeraseLinearRandEffect(filename=oparams_file, **oparams)
    elif model_name == 'TopoILinear':
        my_model = TopoILinear(filename=oparams_file, **oparams)
    elif model_name == 'GyraseLinear':
        my_model = GyraseLinear(filename=oparams_file, **oparams)
    elif model_name == 'TopoISimpleTorque':
        my_model = TopoISimpleTorque(filename=oparams_file, **oparams)
    elif model_name == 'TopoIContinuum':
        my_model = TopoIContinuum(filename=oparams_file, **oparams)
    elif model_name == 'GyraseContinuum':
        my_model = GyraseContinuum(filename=oparams_file, **oparams)
    else:
        raise ValueError('Could not recognise effect model ' + model_name)
    return my_model


# Supercoiling injection of topoisomerases. It injects according the k_cat (injected twist per second), so be careful
# because it can be both positive or negative
def topoisomerase_supercoiling_injection(k_cat, dt):
    position = 0.0
    # Note that k_cat is divided by two on each side because it is assumed that k_cat acts on the local region
    # (both sides)
    twist_left = 0.5 * k_cat * params.w0 * dt
    twist_right = 0.5 * k_cat * params.w0 * dt
    return position, twist_left, twist_right


#  def topoisomerase_lineal_supercoiling_injection(topo, dt):
#    position = 0.0


# ---------------------------------------------------------------------------------------------------------------------
# RNAP FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

# We will use more than once this calculation, so let's store it as function
def uniform_motion(z, dt):
    # Object moves: simple uniform motion

    # Let's get parameters needed
    velocity = z.effect_model.velocity
    gamma = z.effect_model.gamma
    direction = z.direction

    # Calculate change in position.
    position = direction * velocity * dt

    # Injects twist: denatures w = gamma*velocity*dt base-pairs
    twist_left = -direction * gamma * params.w0 * velocity * dt
    twist_right = direction * gamma * params.w0 * velocity * dt
    return position, twist_left, twist_right


# Returns new RNAP position, and the twist it injected on the left and right
def rnap_uniform_motion(z, z_list, dt):
    # Everything 0 for now
    position = 0.0
    twist_left = 0.0
    twist_right = 0.0
    # Get neighbour enzyme
    if z.direction > 0:
        z_n = [e for e in z_list if e.position > z.position][0]  # after - On the right
    if z.direction < 0:
        z_n = [e for e in z_list if e.position < z.position][-1]  # before - On the left
    if z.direction == 0:
        print('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')
        sys.exit()

    # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid clashing
    if z.direction > 0:  # Moving to the right
        if z_n.position - (z.position + position) <= 0:
            return position, twist_left, twist_right
    else:  # Moves to the left
        if z_n.position - (z.position + position) >= 0:
            return position, twist_left, twist_right

    # Nothing is next, so the object moves: simple uniform motion
    position, twist_left, twist_right = uniform_motion(z, dt)

    return position, twist_left, twist_right


# Returns new RNAP position, and the twist it injected on the left and right - It can stall according the Geng
# model of RNAP elongation. In this model, either the RNAP moves with constant velocity or it stalls. It relies on a
# low stretching force params.f_stretching, which here we assume that all molecules interacting exert on the DNA, which
# might not be true... Additionally, if the DNA becomes hyper supercoiled, the RNAP would also stall
def rnap_torque_stall_Geng(z, z_list, dt):
    # For now, nothing happens
    position = 0.0
    twist_left = 0.0
    twist_right = 0.0
    # Get enzymes on the right and left
    z_right = [e for e in z_list if e.position > z.position][0]  # after - On the right
    z_left = [e for e in z_list if e.position < z.position][-1]  # before - On the left
    # Calculate torques and determine if the RNAP will stall
    torque_right = Marko_torque(z.superhelical)  # Torque on the right
    torque_left = Marko_torque(z_left.superhelical)  # Torque on the left
    torque = abs(torque_left - torque_right)
    if torque >= params.stall_torque:  # If torque higher than the stall torque, the RNAP stalls and doesn't move
        return position, twist_left, twist_right

    # Ok, it didn't stall, now we need the neighbours
    if z.direction > 0:
        z_n = z_right  # after - On the right
    if z.direction < 0:
        z_n = z_left  # before - On the left
    if z.direction == 0:
        print('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')
        sys.exit()

    # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid clashing
    if z.direction > 0:  # Moving to the right
        if z_n.position - (z.position + position) <= 0:
            return position, twist_left, twist_right
    else:  # Moves to the left
        if z_n.position - (z.position + position) >= 0:
            return position, twist_left, twist_right

    # It passed all the filters and didn't stall, now, the object moves with simple uniform motion
    # position = Z.position + Z.direction * v0 * dt
    position = z.direction * v0 * dt

    # Injects twist: denatures w=gamma*v0*dt base-pairs
    twist_left = -z.direction * z.k_cat * v0 * dt
    twist_right = z.direction * z.k_cat * v0 * dt

    return position, twist_left, twist_right


# Torque calculated using Marko's elasticity model
def Marko_torque(sigma):
    if np.abs(sigma) <= np.abs(params.sigma_s):
        torque = sigma * params.cs_energy / params.w0
    elif abs(params.sigma_s) < abs(sigma) < abs(params.sigma_p):
        torque = np.sqrt(
            2 * params.p_stiffness * params.g_energy / (1 - params.p_stiffness / params.cs_energy)) / params.w0_nm
    elif abs(sigma) > abs(params.sigma_p):
        torque = sigma * params.p_stiffness / params.w0
    else:
        print('Error in Marko_torque function')
        sys.exit()
    return torque


#  The velocity has the form: v = vmax/ (1+e^{k(T_0 - T_c)} )
#  where vmax = maximum velocity, k = torque parameter, T_0 = Torque acting on enzyme
#  and T_c = cutoff or stalling torque.
#  This function is based on the 2022SevierBioJ paper
def velocity_2022SevierBioJ(z, torque):
    # top = 2.0 * z.effect_model.velocity
    top = z.effect_model.velocity
    exp_arg = z.effect_model.kappa * (torque - z.effect_model.stall_torque)
    exp_arg = np.float128(exp_arg)

    # Define a maximum value for the argument to exp to prevent overflow
    max_exp_arg = 709  # slightly below the overflow threshold for float64
    # Clip the argument to the maximum value
    exp_arg_clipped = np.clip(exp_arg, a_min=None, a_max=max_exp_arg)
    #down = 1.0 + np.exp(exp_arg)
    down = 1.0 + np.exp(exp_arg_clipped)
    velocity = top / down
    return velocity

# ---------------------------------------------------------------------------------------------------------------------
# USEFUL FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
