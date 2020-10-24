"""
Simulations of a flying disc. This file contains the disc object.

Physics are carried out elsewhere.
"""

import numpy as np

from frispy.environment import Environment
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

    def __init__(self, environment=Environment()):
        self.environment = environment

        # Attach an uninitialized trajectory object
        self.initialize_trajectory()

    def initialize_trajectory(self, **kwargs):
        """Set initial conditions for the trajectory."""
        self._trajectory = Trajectory(**kwargs)

    def compute_trajectory(self):
        """Call the differential equation solver to compute
        the trajectory. Return all kinematic variables
        and the timesteps.
        """
        pass

    def set_model(self, **kwargs):
        """Specify the physics model to be used to
        compute the forces and torques.
        """
        self._trajectory.set_model(**kwargs)

    def get_model(self):
        """Return the physics model for the forces
        and torques.
        """
        return self._trajectory.get_model()


if __name__ == "__main__":
    d = Disc()
    print(d)
    print(dir(d))
