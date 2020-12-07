from typing import Dict, Optional

from frispy.equations_of_motion import EOM
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
      eom (EOM, optional): the equations of motion
      initial_conditions (Dict, optional): initial conditions of the disc in
        flight units are in in "mks" and angles are in radians. By defualt,
        the initial conditions will be that the disc is 1 meter off the
        ground (`z=1`), moving at 10 meters/sec in the `x` direction
        (`vx=10`) and is spinning about it's vertical axis at approximately
        10 revolutions per second (approx. 62 rad/sec, or `dgamma=62`).
        All other kinematic variables are set to 0. This configuration
        results in an angle of attack of 0.
    """

    def __init__(
        self,
        eom: EOM = EOM(),
        initial_conditions: Optional[Dict[str, float]] = None,
    ):
        self._eom = eom
        self.set_default_initial_conditions(initial_conditions)
        self.reset_initial_conditions()

    def compute_trajectory(self):
        """Call the differential equation solver to compute
        the trajectory. Return all kinematic variables
        and the timesteps.
        """
        pass

    def reset_initial_conditions(self) -> None:
        """
        Set the initial_conditions of the disc to the default and
        clear the trajectory.
        """
        self.initial_conditions = self.default_initial_conditions
        self.current_coordinates = self.initial_conditions.copy()
        self.current_trajectory = None
        return

    def set_default_initial_conditions(
        self, initial_conditions: Optional[Dict[str, float]]
    ) -> None:
        base_ICs = {
            "x": 0,
            "y": 0,
            "z": 1.0,
            "vx": 10.0,
            "vy": 0,
            "vz": 0,
            "phi": 0,
            "theta": 0,
            "gamma": 0,
            "dphi": 0,
            "dtheta": 0,
            "dgamma": 62.0,
        }
        for i in base_ICs:
            if initial_conditions is not None:
                assert (
                    i in initial_conditions
                ), f"{i} missing from initial conditions"
        self._default_initial_conditions = initial_conditions or base_ICs

    @property
    def default_initial_conditions(self) -> Dict[str, float]:
        return self._default_initial_conditions

    @property
    def environment(self) -> Environment:
        return self._eom.environment

    @property
    def eom(self) -> EOM:
        return self._eom

    @property
    def model(self) -> Model:
        return self._eom.model

    @property
    def trajectory(self) -> Trajectory:
        return self._eom.trajectory
