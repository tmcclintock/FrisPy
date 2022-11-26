"""Physical model for the forces and torques on a disc."""

from dataclasses import dataclass

import numpy as np


@dataclass
class Model:
    """Kinematic model for a disc.

    Holds all of the aerodynamic parameters coupling the kinematic variables
    (velocities, and angular velocities) to the force and torque magnitudes.

    Args:
        PL0: lift force at zero angle of attack (aoa; i.e. the intercept)
        PLa: lift force linear scaling of aoa
        PD0: drag force at zero aoa (i.e. the intercept)
        PDa: drag force linear scaling of aoa
        PTxwx: torque around x-hat linear scaling with x-hat angular velocity
        PTxwz: torque around x-hat linear scaling with z-hat angular velocity (spin)
        PTy0: torque around y-hat intercept
        PTya: torque around y-hat linear scaling with aoa
        PTywy: torque around y-hat linear scaling with y-hat angular velocity
        PTzwz: torque around z-hat linear scaling with z-hat angular velocity (spin)
        alpha_0: angle of attack at minimal drag (i.e. critical aoa)
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
        """Lift force scale factor.

        Linear in the angle of attack (``alpha``).

        Args:
            alpha: angle of attack in radians

        Returns:
            lift force scale factor
        """
        return self.PL0 + self.PLa * alpha

    def C_drag(self, alpha: float) -> float:
        """Drag force scale factor.

        Quadratic in the angle of attack (``alpha``).

        Args:
            alpha: angle of attack in radians

        Returns:
             drag force scale factor
        """
        return self.PD0 + self.PDa * (alpha - self.alpha_0) ** 2

    def C_x(self, wx: float, wz: float) -> float:
        """Torque scale factor around x-hat.

        Linearly in the z angular velocity
        (``w_z``) and the x angular velocity (``w_x``).

        Args:
            wx: x angular velocity in radians per second
            wz: z angular velocity in radians per second

        Returns:
            torque scale factor around x-hat
        """
        return self.PTxwx * wx + self.PTxwz * wz

    def C_y(self, alpha: float, wy: float) -> float:
        """Toruq scalef actor around y-hat.

        Linearly in the y angular velocity
        (``w_y``) and the angle of attack (``alpha``).

        Args:
            alpha: angle of attack in radians
            wy: y angular velocity in radians per second

        Returns:
            torque scale factor around y-hat
        """
        return self.PTy0 + self.PTywy * wy + self.PTya * alpha

    def C_z(self, wz: float) -> float:
        """Torque scale factor around z-hat.

        Linear in the z angular velocity (``w_z``).

        Args:
            wz: z angular velocity in radians per second

        Returns:
            torque scale factor around z-hat
        """
        return self.PTzwz * wz
