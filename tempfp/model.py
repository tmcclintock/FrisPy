"""Physical model for the forces and torques on a disc."""

#import numpy as np

from typing import Dict

class Model():
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

    def set_values(self, **kwargs: Dict[str, float]) -> None:
        """Set the values of the physics model."""

    def func1(self) -> None:
        """A function to write."""
