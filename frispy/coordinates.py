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
    def rotation_matrix(self) -> np.ndarray:
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
