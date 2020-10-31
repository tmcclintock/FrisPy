"""
The ``EOM`` (or "equations of motion)" object is used to actually run the
ODE solver.
"""

from typing import Dict, Union

import numpy as np

from frispy.environment import Environment
from frispy.model import Model
from frispy.trajectory import Trajectory


class EOM:
    """
    ``EOM`` is short for "equations of motion" is used to run the ODE solver
    from `scipy`. It takes in a model for the disc, the trajectory object,
    the environment, and implements the functions for calculating forces
    and torques.
    """

    def __init__(
        self, environment: Environment(), model: Model(), trajectory: Trajectory()
    ):
        self._environment = environment
        self._model = model
        self._trajectory = trajectory

    @property
    def environment(self) -> Environment:
        return self._environment

    @property
    def model(self) -> Model:
        return self._model

    @property
    def trajectory(self) -> Trajectory:
        return self._trajectory

    def compute_forces(
        self,
        phi: float,
        theta: float,
        velocity: np.ndarray,
        ang_velocity: np.ndarray,
    ) -> Dict[str, Union[float, np.ndarray, Dict[str, np.ndarray]]]:
        """
        Compute the lift, drag, and gravitational forces on the disc.

        Args:
        TODO
        """
        res = self.trajectory.calculate_intermediate_quantities(
            phi, theta, velocity, ang_velocity
        )
        aoa = res["angle_of_attack"]
        vhat = velocity / np.linalg.norm(velocity)
        force_amplitude = (
            0.5
            * self.environment["air_density"]
            * (velocity @ velocity)
            * self.environment["mass"]
        )
        # Compute the lift and drag forces
        res["F_lift"] = (
            self.model.C_lift(aoa)
            * force_amplitude
            * np.cross(vhat, res["unit_vectors"]["yhat"])
        )
        res["F_drag"] = self.model.C_drag(aoa) * force_amplitude * (-vhat)
        # Compute gravitational force
        res["F_grav"] = (
            self.environment["mass"]
            * self.environment["g"]
            * self.environment["grav_vector"]
        )
        return res

    def compute_torques(
        self,
        velocity: np.ndarray,
        ang_velocity: np.ndarray,
        res: Dict[str, Union[float, np.ndarray, Dict[str, np.ndarray]]],
    ) -> Dict[str, Union[float, np.ndarray, Dict[str, np.ndarray]]]:
        """
        Compute the torque around each principle axis.

        Args:
        TODO
        """
        aoa = res["angle_of_attack"]
        torque_amplitude = (
            0.5
            * self.environment["air_density"]
            * (velocity @ velocity)
            * self.environment["diameter"]
            * self.environment["area"]
        )
        wx, wy, wz = res["w"]
        # Compute component torques. Note that "x" and "y" are computed
        # in the lab frame
        res["T_x_lab"] = (
            self.model.C_x(wx, wz)
            * torque_amplitude
            * res["unit_vectors"]["xhat"]
        )
        res["T_y_lab"] = (
            self.model.C_y(aoa, wy)
            * torque_amplitude
            * res["unit_vectors"]["yhat"]
        )
        # Computed in the disc frame
        res["T_z"] = self.model.C_z(wz) * torque_amplitude * np.array([0, 0, 1])
        # Rotate into the disc frame
        res["T_x"] = res["rotation_matrix"] @ res["T_x_lab"]
        res["T_y"] = res["rotation_matrix"] @ res["T_y_lab"]
        res["T"] = res["T_x"] + res["T_y"] + res["T_z"]
        return res
