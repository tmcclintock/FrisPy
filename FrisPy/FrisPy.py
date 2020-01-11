import numpy as np

from .trajectory import Trajectory
#from .disc import * #not PEP8 compliant

class Disc(object):
    """Flying spinning disc object. The disc object contains only physical parameters
    of the disc and environment that it exists (e.g. gravitational acceleration and
    air density). Note that the default area, mass, and inertial moments are for
    Discraft Ultrastars (175 grams or 0.175 kg).

    All masses are kg, lengths are meters (m) and times are seconds (s). That is,
    these files all use `mks` units.

    air_density (float): default is 1.225 kg/m^3
    area (float): default is 0.057 m^2
    diameter (float): default is 0.27 m; (computed from default area)
    g (float): default is 9.81 m/s^2; gravitational acceleration on Earth
    grav_vector (numpy.ndarray): default is [0,0,-1]; is a unit vector
    I_zz (float): default is 0.001219 kg*m^2; z-axis moment of inertia
    I_xx (float): default is 0.175 kg*m^2; x and y-axis moments of inertia
        (i.e. is the same as I_yy and the cross components I_xy)
    mass (float): defualt is 0.175 kg

    """
    def __init__(self, **kwargs):
        # all those keys will be initialized as class attributes
        allowed_keys = ['air_density', 'area', 'diameter', 'g',
                        'grav_vector', 'I_zz', 'I_xx', 'mass']
        default_values = [1.225, 0.057, 9.81,
                          np.array([0,0,-1]), 0.002352, 0.001219, 0.175]
        default_values.insert(2, 2*(default_values[1]/np.pi)**0.5)
        # initialize all allowed keys to defaults
        self.__dict__.update((key, value) for key, value in zip(allowed_keys, default_values))
        # and update the given keys by their given values
        self.__dict__.update((key, value) for key, value in kwargs.items() if key in allowed_keys)

        #Attach an uninitialized trajectory object
        self._trajectory = Trajectory()

    def initialize_trajectory(self):
        pass

    def compute_trajectory(self):
        pass

    def set_model(self):
        pass
        
if __name__ == "__main__":
    d = Disc()
    print(d)
    print(dir(d))
    print(d.mass)
