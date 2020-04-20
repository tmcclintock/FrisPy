"""An object for holding an interface to the differential equation solver
for the disc trajectory.
"""

#import numpy as np

from .model import Model

class Trajectory():
    """Trajectory of a disc object. This object contains the initial
    values, and computed values for the kinematic variables of
    the disc including the position, velocities, angles, and
    angular velocities. Units are all given in `FrisPy.py` in
    the documentation for the `Disc` object.

    Args:
        x (float): horizontal position; default is 0 m
        y (float): horizontal position; default is 0 m
        z (float): vertical position; default is 1 m
        vx (float): x-velocity; default is 10 m/s
        vy (float): y-velocity; default is 0 m/s
        vz (float): z-velocity; default is 0 m/s
        phi (float): 1st Euler angle (pitch); default is 0 rad
        theta (float): 2nd Euler angle (roll); default is 0 rad
        gamma (float): 3rd Euler angle (spin); default is 0 rad
        phidot (float): phi angular velocity; default is 0 rad/s
        thetadot (float): theta angular velocity; default is 0 rad/s
        gammadot (float): gamma angular velocity; default is 50 rad/s

    """
    def __init__(self, **kwargs):
        allowed_keys = ['x', 'y', 'z', 'vx', 'vy', 'vz',
                        'phi', 'theta', 'gamma',
                        'phidot', 'thetadot', 'gammadot']
        default_values = [0, 0, 1, 10, 0, 0,
                          0, 0, 0, 0, 0, 50]
        # initialize all allowed keys to defaults
        self.__dict__.update((key, value) for key, value in
                             zip(allowed_keys, default_values))
        # and update the given keys by their given values
        self.__dict__.update((key, value) for key, value in
                             kwargs.items() if key in allowed_keys)

        self.initial_values = {"initial_"+key: value for key, value in
                               self.__dict__.items() if key in allowed_keys}

        self._model = Model()

    def set_model(self, **kwargs):
        """Update the values of parameters in the
        coefficient model.

        Args:
            kwargs (dict): key-value pairs of coefficients in the model

        Returns:
            None
        """
        self._model.set_values(**kwargs)

    def get_model(self):
        """Return the model object."""
        return self._model
