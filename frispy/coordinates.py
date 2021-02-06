"""
Coordinates of the disc. This includes all kinematic variables
as well as positions (i.e. euclidean position, velocities,
angular coordinates, angular velocities).
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Coordinates:
    __slots__ = [
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "_phi",
        "_theta",
        "gamma",
        "dphi",
        "dtheta",
        "dgamma",
        "_sin_phi",
        "_cos_phi",
        "_sin_theta",
        "_cos_theta",
    ]
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    _phi: float
    _theta: float
    gamma: float
    dphi: float
    dtheta: float
    dgamma: float

    def __init__(
        self,
        x: float,
        y: float,
        z: float,
        vx: float,
        vy: float,
        vz: float,
        phi: float,
        theta: float,
        gamma: float,
        dphi: float,
        dtheta: float,
        dgamma: float,
    ):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.phi = phi
        self.theta = theta
        self.gamma = gamma
        self.dphi = dphi
        self.dtheta = dtheta
        self.dgamma = dgamma

    @property
    def phi(self):
        return self._phi

    @property
    def sin_phi(self):
        return self._sin_phi

    @property
    def cos_phi(self):
        return self._cos_phi

    @phi.setter
    def phi(self, value):
        self._sin_phi = np.sin(value)
        self._cos_phi = np.cos(value)
        self._phi = value

    @property
    def theta(self):
        return self._theta

    @property
    def sin_theta(self):
        return self._sin_theta

    @property
    def cos_theta(self):
        return self._cos_theta

    @theta.setter
    def theta(self, value):
        self._sin_theta = np.sin(value)
        self._cos_theta = np.cos(value)
        self._theta = value

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
        sp = self.sin_phi
        cp = self.cos_phi
        st = self.sin_theta
        ct = self.cos_theta
        return np.array([[ct, sp * st, -st * cp], [0, cp, sp], [st, -sp * ct, cp * ct]])
