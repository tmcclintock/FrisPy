import numpy as np

#from .disc import * #not PEP8 compliant

class Disc(object):
    def __init__(self, **kwargs):
        # all those keys will be initialized as class attributes
        allowed_keys = ['air_density', 'area', 'diameter', 'g',
                        'grav_vector', 'I_zz', 'I_xx', 'mass']
        default_values = [1.225, 0.057, 2*(0.057/np.pi)**0.5, 9.81,
                          np.array([0,0,-1]), 0.002352, 0.001219, 0.175]
        # initialize all allowed keys to false
        self.__dict__.update((key, value) for key, value in zip(allowed_keys, default_values))
        # and update the given keys by their given values
        self.__dict__.update((key, value) for key, value in kwargs.items() if key in allowed_keys)


if __name__ == "__main__":
    d = Disc()
    print(d)
    print(dir(d))
