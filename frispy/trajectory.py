"""An object for holding an interface to the differential equation solver
for the disc trajectory.
"""

from numbers import Number
from typing import Dict

import numpy as np


class Trajectory:
    """
    Class for computing the disc flight trajectory. Takes initial values
    and interfaces with an ODE solver.

    Units are meters [m] for length, kilograms [kg] for mass, seconds [s]
    for time, and radians [rad] for angles.

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
        # A default flight configuration
        self._initial_conditions: Dict[str, float] = {
            "x": 0,
            "y": 0,
            "z": 1,
            "vx": 10,
            "vy": 0,
            "vz": 0,
            "phi": 0,
            "theta": 0,
            "gamma": 0,
            "phidot": 0,
            "thetadot": 0,
            "gammadot": 50,
        }
        self._coord_order = [
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
            "phi",
            "theta",
            "gamma",
            "phidot",
            "thetadot",
            "gammadot",
        ]

        # set arguments to initial conditions
        for k, v in kwargs.items():
            assert (
                k in self._initial_coordinates
            ), f"invalid initial condition name {k}"
            assert isinstance(v, Number), f"invalid type for {v}, {type(v)}"
            self._initial_conditions[k] = v

    @property
    def initial_conditions(self) -> Dict[str, float]:
        return self._initial_conditions

    @property
    def initial_conditions_array(self) -> np.ndarray:
        return np.array([self.initial_conditions[k] for k in self._coord_order])

    # @property
    # def velocity(self) -> np.ndarray:
    #     return np.array([
    # ])

    @staticmethod
    def rotation_matrix(phi: float, theta: float) -> np.ndarray:
        """
        Compute the (partial) rotation matrix that transforms from the
        lab frame to the disc frame. Note that because of azimuthal
        symmetry, the azimuthal angle (`gamma`) is not used.
        """
        sp = np.sin(phi)
        cp = np.cos(phi)
        st = np.sin(theta)
        ct = np.cos(theta)
        return np.array(
            [[ct, sp * st, -st * cp], [0, cp, sp], [st, -sp * ct, cp * ct]]
        )
