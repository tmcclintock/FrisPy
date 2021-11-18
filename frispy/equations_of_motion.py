from numbers import Number
from typing import Dict, Union

import numpy as np

from frispy.environment import Environment
from frispy.model import Model


class EOM:
    """
    ``EOM`` is short for "equations of motion" is used to run the ODE solver
    from `scipy`. It takes in a model for the disc, the trajectory object,
    the environment, and implements the functions for calculating forces
    and torques.
    """

    def __init__(
        self,
        area: Number,
        I_xx: Number,
        I_zz: Number,
        mass: Number,
        environment: Environment = Environment(),
        model: Model = Model(),
    ):
        self.area = area
        self.diameter = 2 * np.sqrt(self.area / np.pi)
        self.I_xx = I_xx
        self.I_zz = I_zz
        self.mass = mass
        self.environment = environment
        self.model = model
        # Pre-compute some values to optimize the ODEs
        self.force_per_v2 = (
            0.5 * self.environment.air_density * self.area
        )  # N / (m/s)^2
        self.torque_per_v2 = self.force_per_v2 * self.diameter  # N * m / (m/s)^2
        self.F_grav = (
            self.mass * self.environment.g * self.environment.grav_unit_vector
        )  # looks like [0, 0, -m*g]
        self.z_hat = np.array([0, 0, 1])

    @classmethod
    def rotation_matrix_from_phi_theta(
        cls, phi: float, theta: float
    ) -> np.ndarray:
        sp, cp = np.sin(phi), np.cos(phi)
        st, ct = np.sin(theta), np.cos(theta)
        return cls.rotation_matrix(sp, cp, st, ct)

    @staticmethod
    def rotation_matrix(sp: float, cp: float, st: float, ct) -> np.ndarray:
        """
        Compute the (partial) rotation matrix that transforms from the
        lab frame to the disc frame. Note that because of azimuthal
        symmetry, the azimuthal angle (`gamma`) is not used.

        This matrix (R) can be used to transform a vector from the lab frame (L)
        into the disk frame (D), i.e.: r_D = R dot r_L.

        The ``z_hat`` unit vector in the disk frame (D) will always be pointing
        perpendicular up from the top face of the disk.
        """
        return np.array(
            [[ct, sp * st, -st * cp], [0, cp, sp], [st, -sp * ct, cp * ct]]
        )

    @classmethod
    def compute_angle_of_attack(
        cls,
        phi: float,
        theta: float,
        velocity: np.ndarray,
        return_all_variables: bool = False,
    ):
        # Rotation matrix
        sp, cp = np.sin(phi), np.cos(phi)
        st, ct = np.sin(theta), np.cos(theta)
        rotation_matrix = cls.rotation_matrix(sp, cp, st, ct)
        # Unit vectors
        zhat = rotation_matrix[2]
        v_dot_zhat = velocity @ zhat
        v_in_plane = velocity - zhat * v_dot_zhat
        angle_of_attack = -np.arctan(v_dot_zhat / np.linalg.norm(v_in_plane))
        if return_all_variables:
            return (
                angle_of_attack,
                sp,
                cp,
                st,
                ct,
                rotation_matrix,
                v_dot_zhat,
                v_in_plane,
            )
        else:
            return angle_of_attack

    def geometric_quantities(
        self,
        phi: float,
        theta: float,
        velocity: np.ndarray,
        angular_velocity: np.ndarray,
    ) -> Dict[str, Union[float, np.ndarray, Dict[str, np.ndarray]]]:
        """
        Compute intermediate quantities on the way to computing the time
        derivatives of the kinematic variables.
        """
        (
            angle_of_attack,
            _,
            _,
            st,
            ct,
            rotation_matrix,
            _,
            v_in_plane,
        ) = self.compute_angle_of_attack(
            phi, theta, velocity, return_all_variables=True
        )
        zhat: np.ndarray = rotation_matrix[2]
        xhat: np.ndarray = v_in_plane / np.linalg.norm(v_in_plane)
        yhat = np.cross(zhat, xhat)
        # Angle of attack
        # Disc angular velocities
        angular_velocity = angular_velocity
        w_prime = np.array(
            [
                angular_velocity[0] * ct,
                angular_velocity[1],
                angular_velocity[0] * st + angular_velocity[2],
            ]
        )
        # Angular velocity in lab coordinates
        w_lab = w_prime @ rotation_matrix
        # Angular velocity components along the unit vectors in the lab frame
        w = np.array([xhat, yhat, zhat]) @ w_lab
        return {
            "unit_vectors": {"xhat": xhat, "yhat": yhat, "zhat": zhat},
            "angle_of_attack": angle_of_attack,
            "rotation_matrix": rotation_matrix,
            "w_prime": w_prime,
            "w_lab": w_lab,
            "w": w,
        }

    def compute_forces(
        self,
        phi: float,
        theta: float,
        velocity: np.ndarray,
        ang_velocity: np.ndarray,
    ) -> Dict[str, Union[float, np.ndarray, Dict[str, np.ndarray]]]:
        """
        Compute the lift, drag, and gravitational forces on the disc.
        """
        res = self.geometric_quantities(phi, theta, velocity, ang_velocity)
        aoa = res["angle_of_attack"]
        vhat = velocity / np.linalg.norm(velocity)
        force_amplitude = self.force_per_v2 * (velocity @ velocity)
        # Compute the lift and drag forces
        res["F_lift"] = (
            self.model.C_lift(aoa)
            * force_amplitude
            * np.cross(vhat, res["unit_vectors"]["yhat"])
        )
        res["F_drag"] = self.model.C_drag(aoa) * force_amplitude * (-vhat)
        # Compute gravitational force
        res["F_grav"] = self.F_grav
        res["F_total"] = res["F_lift"] + res["F_drag"] + res["F_grav"]
        res["Acc"] = res["F_total"] / self.mass
        return res

    def compute_torques(
        self,
        velocity: np.ndarray,
        res: Dict[str, Union[float, np.ndarray, Dict[str, np.ndarray]]],
    ) -> Dict[str, Union[float, np.ndarray, Dict[str, np.ndarray]]]:
        """
        Compute the torque around each principle axis.
        """
        aoa = res["angle_of_attack"]
        res["torque_amplitude"] = self.torque_per_v2 * (velocity @ velocity)
        wx, wy, wz = res["w"]
        # Compute component torques. Note that "x" and "y" are computed
        # in the lab frame
        res["T_x_lab"] = (
            self.model.C_x(wx, wz)
            * res["torque_amplitude"]
            * res["unit_vectors"]["xhat"]
        )
        res["T_y_lab"] = (
            self.model.C_y(aoa, wy)
            * res["torque_amplitude"]
            * res["unit_vectors"]["yhat"]
        )
        # Computed in the disc frame
        res["T_z"] = self.model.C_z(wz) * res["torque_amplitude"] * self.z_hat
        # Rotate into the disc frame
        res["T_x"] = res["rotation_matrix"] @ res["T_x_lab"]
        res["T_y"] = res["rotation_matrix"] @ res["T_y_lab"]
        res["T"] = res["T_x"] + res["T_y"] + res["T_z"]
        return res

    def compute_derivatives(
        self, time: float, coordinates: np.ndarray
    ) -> np.ndarray:
        """
        Right hand side of the ordinary differential equations. This is
        supplied to :meth:`scipy.integrate.solve_ivp`. See `this page
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_
        for more information about its `fun` argument.

        .. todo::

           Implement the disc hitting the ground as a (Callable) scipy
           event object.

        Args:
          time (float): instantanious time of the system
          coordinates (np.ndarray): kinematic variables of the disc

        Returns:
          derivatives of all coordinates
        """
        # If the disk hit the ground, then stop calculations
        if coordinates[2] <= 0:
            return coordinates * 0

        velocity = coordinates[3:6]
        ang_velocity = coordinates[9:12]
        result = self.compute_forces(
            coordinates[6], coordinates[7], velocity, ang_velocity
        )
        result = self.compute_torques(velocity, result)
        return np.concatenate(
            (velocity, result["Acc"], ang_velocity, result["T"])
        )
