"""
Simulations of a flying disc. This file contains the disc object.

Physics are carried out elsewhere.
"""

import numpy as np

from frispy.environment import Environment
from frispy.model import Model
from frispy.trajectory import Trajectory


class Disc:
    """Flying spinning disc object. The disc object contains only physical
    parameters of the disc and environment that it exists (e.g. gravitational
    acceleration and air density). Note that the default area, mass, and
    inertial moments are for Discraft Ultrastars (175 grams or 0.175 kg).

    All masses are kg, lengths are meters (m) and times are seconds (s). That
    is, these files all use `mks` units. Angular units use radians (rad), and
    angular velocities are in rad/s.

    Args:
        environment (Environment): experiment environment that contains e.g.
            air density, disc mass, etc. See the ``Environment``. Default is
    """

    def __init__(
            self,
            environment=Environment(),
            model=Model(),
            trajectory=Trajectory(),
    ):
        self._environment = environment
        self._model = model
        self._trajectory = trajectory

    def compute_trajectory(self):
        """Call the differential equation solver to compute
        the trajectory. Return all kinematic variables
        and the timesteps.
        """
        pass

    @property
    def environment(self) -> Environment:
        return self._environment

    @property
    def model(self) -> Model:
        return self._model

    @property
    def trajectory(self) -> Trajectory:
        return self._trajectory
