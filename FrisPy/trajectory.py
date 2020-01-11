import numpy as np

class Trajectory(object):
    """Trajectory of a disc object. This object contains the initial
    values, and computed values for the kinematic variables of
    the disc including the position, velocities, angles, and
    angular velocities. Units are all given in `FrisPy.py` in
    the documentation for the `Disc` object.

    

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
