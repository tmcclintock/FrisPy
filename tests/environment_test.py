"""
Tests of the ``Environment`` object.
"""
from unittest import TestCase

import numpy as np

from frispy import Environment


class TestEnvironment(TestCase):
    def test_smoke(self):
        e = Environment()
        assert e is not None

    def test_default_values(self):
        e = Environment()
        d = {
            "air_density": 1.225,
            "area": 0.57,
            "g": 9.81,
            "grav_unit_vector": np.array([0.0, 0.0, -1.0]),
            "I_zz": 0.002352,
            "I_xx": 0.001219,
            "mass": 0.175,
        }
        for k, v in d.items():
            np.testing.assert_equal(getattr(e, k), v, k)
        assert (e.diameter / 2) ** 2 * np.pi == e.area
