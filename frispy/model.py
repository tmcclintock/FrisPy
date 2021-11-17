"""
Physical model for the forces and torques on a disc.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Model:
    """
    Coefficient model for a disc. Holds all of the aerodynamic
    parameters coupling the kinematic variables (spins and angles)
    to the force magnitudes.
    """

    PL0: float = 0.33
    PLa: float = 1.9
    PD0: float = 0.18
    PDa: float = 0.69
    PTxwx: float = -0.013
    PTxwz: float = -0.0017
    PTy0: float = -0.082
    PTya: float = 0.43
    PTywy: float = -0.014
    PTzwz: float = -0.000034
    alpha_0: float = 4 * np.pi / 180

    #####################################################################
    # Below are functions connecting physical variables to force/torque #
    # scaling factors (the `C`s)                                        #
    #####################################################################

    def C_lift(self, alpha: float) -> float:
        """
        Lift force scale factor. Linear in the angle of attack (`alpha`).

        Args:
            alpha (float): angle of attack in radians

        Returns:
            (float) lift force scale factor
        """
        return self.PL0 + self.PLa * alpha

    def C_drag(self, alpha: float) -> float:
        """
        Drag force scale factor. Quadratic in the angle of attack (`alpha`).

        Args:
            alpha (float): angle of attack in radians

        Returns:
            (float) drag force scale factor
        """
        return self.PD0 + self.PDa * (alpha - self.alpha_0) ** 2

    def C_x(self, wx: float, wz: float) -> float:
        """
        'x'-torque scale factor. Linearly additive in the 'z' angular velocity
        (`w_z`) and the 'x' angular velocity (`w_x`).

        Args:
            wx (float): 'x' angular velocity in radians per second
            wz (float): 'z' angular velocity in radians per second

        Returns:
            (float) 'x'-torque scale factor
        """
        return self.PTxwx * wx + self.PTxwz * wz

    def C_y(self, alpha: float, wy: float) -> float:
        """
        'y'-torque scale factor. Linearly additive in the 'y' angular velocity
        (`w_y`) and the angle of attack (`alpha`).

        Args:
            alpha (float): angle of attack in radians
            wy (float): 'y' angular velocity in radians per second

        Returns:
            (float) 'y'-torque scale factor
        """
        return self.PTy0 + self.PTywy * wy + self.PTya * alpha

    def C_z(self, wz: float) -> float:
        """
        'z'-torque scale factor. Linear in the 'z' angular velocity
        (`w_z`).

        Args:
            wz (float): 'z' angular velocity in radians per second

        Returns:
            (float) 'z'-torque scale factor
        """
        return self.PTzwz * wz
