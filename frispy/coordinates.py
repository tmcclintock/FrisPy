"""
Coordinates of the disc. This includes all kinematic variables
as well as positions (i.e. euclidean position, velocities,
angular coordinates, angular velocities).
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Coordinates:
    __slots__ = ["x", "y", "z", "vx", "vy", "vz", "phi", "theta", "gamma", "dphi", "dtheta", "dgamma"]
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    phi: float
    theta: float
    gamma: float
    dphi: float
    dtheta: float
    dgamma: float

    @property
    def position(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

    @property
    def velocity(self) -> np.ndarray:
        return np.array([self.vx, self.vy, self.vz])

    @property
    def angles(self) -> np.ndarray:
        return np.array([self.phi, self.theta, self.gamma])

    @property
    def angular_velocity(self) -> np.ndarray:
        return np.array([self.dphi, self.dtheta, self.dgamma])

    @property
    def rotation_matrix(self) -> np.ndarray:
        """
        Compute the (partial) rotation matrix that transforms from the
        lab frame to the disc frame. Note that because of azimuthal
        symmetry, the azimuthal angle (`gamma`) is not used.
        """
        sp = np.sin(self.phi)
        cp = np.cos(self.phi)
        st = np.sin(self.theta)
        ct = np.cos(self.theta)
        return np.array(
            [[ct, sp * st, -st * cp], [0, cp, sp], [st, -sp * ct, cp * ct]]
        )
