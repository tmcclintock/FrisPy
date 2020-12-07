"""
Physical model for the forces and torques on a disc.
"""

from typing import Dict

import numpy as np


class Model:
    """
    Coefficient model for a disc. Holds all of the aerodynamic
    parameters coupling the kinematic variables (spins and angles)
    to the force magnitudes.
    """

    def __init__(self, **kwargs):
        self._coefficients: Dict[str, float] = {
            "PL0": 0.33,
            "PLa": 1.9,
            "PD0": 0.18,
            "PDa": 0.69,
            "PTxwx": 0.43,
            "PTxwz": -1.4e-2,
            "PTy0": -8.2e-2,
            "PTya": -1.2e-2,
            "PTywy": -1.7e-3,
            "PTzwz": -3.4e-5,
            "alpha_0": 4 * np.pi / 180,
        }
        for k, v in kwargs.items():
            assert k in self.coefficients, f"invalid coefficient name {k}"
            self.coefficients[k] = v

    def set_value(self, name: str, value: float) -> None:
        """
        Set the value of a coefficient.

        Args:
            name (str): name of the coefficient
            value (float): value of the coefficient
        """
        assert name in self.coefficients, f"invalid coefficient name {name}"
        self._coefficients[name] = value

    def set_values(self, coefs: Dict[str, float]) -> None:
        """
        Set the values of the coefficients.

        Args:
            coefs (Dict[str, float]): key-value pairs of coeffient names
                ane their values
        """
        for k, v in coefs.items():
            self.set_value(k, v)

    def get_value(self, name: str) -> float:
        """
        Obtain the value of the coefficient.

        Args:
            name (str): name of the coefficient

        Returns:
            value of the coefficient with the name `name`
        """
        assert name in self.coefficients, f"invalid coefficient name {name}"
        return self.coefficients[name]

    @property
    def coefficients(self) -> Dict[str, float]:
        return self._coefficients

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
        PL0 = self.get_value("PL0")
        PLa = self.get_value("PLa")
        return PL0 + PLa * alpha

    def C_drag(self, alpha: float) -> float:
        """
        Drag force scale factor. Quadratic in the angle of attack (`alpha`).

        Args:
            alpha (float): angle of attack in radians

        Returns:
            (float) drag force scale factor
        """
        PD0 = self.get_value("PD0")
        PDa = self.get_value("PDa")
        alpha_0 = self.get_value("alpha_0")
        return PD0 + PDa * (alpha - alpha_0) ** 2

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
        PTxwx = self.get_value("PTxwx")
        PTxwz = self.get_value("PTxwz")
        return PTxwx * wx + PTxwz * wz

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
        PTy0 = self.get_value("PTy0")
        PTywy = self.get_value("PTywy")
        PTya = self.get_value("PTya")
        return PTy0 + PTywy * wy + PTya * alpha

    def C_z(self, wz: float) -> float:
        """
        'z'-torque scale factor. Linear in the 'z' angular velocity
        (`w_z`).

        Args:
            wz (float): 'z' angular velocity in radians per second

        Returns:
            (float) 'z'-torque scale factor
        """
        PTzwz = self.get_value("PTzwz")
        return PTzwz * wz
