"""
The ``Environment`` object.
"""

from typing import NamedTuple

import numpy as np


class Environment(NamedTuple):
    """
    The environment in which the disc is flying in. This object contains
    information on the magnitude and direction of gravity, properties of the
    wind, and also intrinsic properties of the disc such as its area and
    mass.

    Args:
        air_density (float): default is 1.225 kg/m^3
        area (float): default is 0.057 m^2
        g (float): default is 9.81 m/s^2; gravitational acceleration on Earth
        I_zz (float): default is 0.002352 kg*m^2; z-axis moment of inertia
        I_xx (float): default is 0.001219 kg*m^2; x and y-axis moments of inertia
            (i.e. is the same as I_yy and the cross components I_xy)
        mass (float): defualt is 0.175 kg
    """

    air_density: float = 1.225  # kg/m^3
    area: float = 0.058556  # m^2
    g: float = 9.81  # m/s^2
    I_zz: float = 0.002352  # kg*m^2
    I_xx: float = 0.001219  # kg*m^2
    mass: float = 0.175  # kg

    @property
    def grav_unit_vector(self) -> np.ndarray:
        """Gravitational direction."""
        return np.array([0, 0, -1])

    @property
    def diameter(self) -> float:
        """Disc diameter."""
        return 2 * (self.area / np.pi) ** 0.5
