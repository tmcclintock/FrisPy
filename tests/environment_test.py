"""
Tests of the ``Environment`` object.
"""
import numpy as np

from frispy import Environment


def test_smoke():
    e = Environment()
    assert e is not None


def test_default_values():
    e = Environment()
    d = {
        "air_density": 1.225,
        "g": 9.81,
        "grav_unit_vector": np.array([0.0, 0.0, -1.0]),
    }
    for k, v in d.items():
        np.testing.assert_equal(getattr(e, k), v, k)
