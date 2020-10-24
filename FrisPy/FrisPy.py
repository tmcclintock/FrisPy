"""
Simulations of a flying disc. This file contains the disc object.

Physics are carried out elsewhere.
"""

import numpy as np

from trajectory import Trajectory

# from .disc import * #not PEP8 compliant


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
