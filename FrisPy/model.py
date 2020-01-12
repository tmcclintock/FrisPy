import numpy as np

class Model(object):
    """Coefficient model for a disc. Holds all of the aerodynamic
    parameters coupling the kinematic variables (spins and angles)
    to the force magnitudes.

    """
    def __init__(self, **kwargs):
        allowed_keys = ["a"]
        default_values = [0]
        # initialize all allowed keys to defaults
        self.__dict__.update((key, value) for key, value in
                             zip(allowed_keys, default_values))
        # and update the given keys by their given values
        self.__dict__.update((key, value) for key, value in
                             kwargs.items() if key in allowed_keys)

