"""
The ``Environment`` object.
"""

import numpy as np


class Environment(dict):
    """
    The environment in which the disc is flying in. This object contains
    information on the magnitude and direction of gravity, properties of the
    wind, and also intrinsic properties of the disc such as its area and
    mass.

    Args:
        air_density (float): default is 1.225 kg/m^3
        area (float): default is 0.057 m^2
        g (float): default is 9.81 m/s^2; gravitational acceleration on Earth
        grav_vector (List[float]): default is [0,0,-1]; is a unit vector
        I_zz (float): default is 0.001219 kg*m^2; z-axis moment of inertia
        I_xx (float): default is 0.175 kg*m^2; x and y-axis moments of inertia
            (i.e. is the same as I_yy and the cross components I_xy)
        mass (float): defualt is 0.175 kg
    """

    def __init__(self, **kwargs):
        super().__init__(
            {
                "air_density": 1.225,
                "area": 0.57,
                "g": 9.81,
                "grav_vector": np.array([0.0, 0.0, -1.0]),
                "I_zz": 0.002352,
                "I_xx": 0.001219,
                "mass": 0.175,
            }
        )
        # Insert user-defined quantities
        for k, v in kwargs.items():
            self[k] = v

        # Set diameter by hand from area
        if "diameter" not in self:
            self["diameter"] = 2 * (self["area"] / np.pi) ** 0.5
        np.testing.assert_almost_equal(
            self["area"] / (self["diameter"] / 2) ** 2,
            np.pi,
            err_msg=f"area {self['area']} and diameter {self['diameter']} "
            + "not geometrically consistent",
        )
