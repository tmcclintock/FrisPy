"""
Physical model for the forces and torques on a disc.
"""

from typing import Dict

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
            "PTxwx":0.43,
            "PTxwz":-1.4e-2,
            "PTy0":-8.2e-2,
            "PTya":-1.2e-2,
            "PTywy":-1.7e-3,
            "PTzwz":-3.4e-5,
        }
        for k, v in kwargs.items():
            assert k in self._coefficients, f"invalid coefficient name {k}"
            self._coefficients[k] = v

    def set_values(self, coefs: Dict[str, float] ) -> None:
        """
        Set the values of the coefficients.

        Args:
            coefs (Dict[str, float]): key-value pairs of coeffient names
                ane their values
        """
        for k, v in coefs.items():
            assert k in self._coefficients, f"invalid coefficient name {k}"
            self._coefficients[k] = v

    def get_value(self, name: str) -> float:
        """
        Obtain the value of the coefficient.

        Args:
            name (str): name of the coefficient

        Returns:
            value of the coefficient with the name `name`
        """
        assert name in self._coefficients, f"invalid coefficient name {k}"
        return self._coefficients[name]

    def func1(self) -> None:
        """A function to write."""
