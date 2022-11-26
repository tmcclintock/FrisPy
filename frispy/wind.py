"""
The ``Wind`` class handles the wind, which is a real-valued time-dependent
vector field that influences the flight of the disc.
"""

from abc import ABC, abstractmethod
from typing import List, Optional, Union

import numpy as np


class Wind(ABC):
    """
    Abstract class to handle different types of wind. These can include
    steady, laminar flow winds or swirling winds. Winds can have a time
    dependence to mimic "gusts".
    """

    @abstractmethod
    def wind_vector(
        self,
        t: Optional[Union[float, int, np.ndarray]],
        position: Optional[Union[List, np.ndarray]],
    ) -> np.ndarray:
        """
        Obtain a length 3 vector of the wind at time `t`.
        """
        return np.array([0, 0, 0])


class NoWind(Wind):
    """
    No wind.
    """

    def wind_vector(self, *args) -> np.ndarray:
        """
        All components are zero.
        """
        return np.array([0, 0, 0])


class ConstantWind(Wind):
    """
    The wind is uniform in position and constant in time.
    """

    def __init__(
        self,
        wind_vector: Optional[np.ndarray] = None,
    ):
        """Constructor."""
        super().__init__()
        self._wind_vector = wind_vector or np.zeros(3)

    def wind_vector(self, *args) -> np.ndarray:
        """Return ``self._wind_vector``."""
        return self._wind_vector
