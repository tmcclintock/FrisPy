"""
The ``Wind`` class handles the wind, which is a real-valued time-dependent
vector field that influences the flight of the disc.
"""

from typing import List, Optional, Union

from abc import ABC, abstractmethod

import numpy as np


class Wind(ABC):
    """
    Abstract class to handle different types of wind. These can include
    steady, laminar flow winds or swirling winds. Winds can have a time
    dependence to mimic "gusts".
    """

    @abstractmethod
    def get_wind_vector(
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

    def get_wind(self, *args) -> np.ndarray:
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
        vx: Optional[float] = None,
        vy: Optional[float] = None,
        vz: Optional[float] = None,
        wind_vector: Optional[np.ndarray] = None,
    ):
        super().__init__()
        assert any(
            (vx, vy, vz, wind_vector)
        ), "must provide at least one wind component"
        if wind_vector is not None:
            assert not any(
                (vx, vy, vz)
            ), "cannot provide `wind_vector` and individual components"
            self._wind_vector = wind_vector
        else:
            self._wind_vector = np.zeros([vx or 0, vy or 0, vz or 0])

    def get_wind_vector(self, *args) -> np.ndarray:
        return self._wind_vector
