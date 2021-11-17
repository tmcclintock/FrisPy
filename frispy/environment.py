from dataclasses import dataclass

import numpy as np


@dataclass
class Environment:
    """
    The environment in which the disc is flying in. This object contains
    information on the magnitude and direction of gravity, properties of the
    wind, and also intrinsic properties of the disc such as its area and
    mass.

    Args:
        air_density (float): default is 1.225 kg/m^3
        g (float): default is 9.81 m/s^2; gravitational acceleration on Earth
    """

    air_density: float = 1.225  # kg/m^3
    g: float = 9.81  # m/s^2

    @property
    def grav_unit_vector(self) -> np.ndarray:
        """Gravitational direction."""
        return np.array([0, 0, -1])
