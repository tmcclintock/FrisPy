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
        "area": 0.058556,
        "g": 9.81,
        "grav_unit_vector": np.array([0.0, 0.0, -1.0]),
        "I_zz": 0.002352,
        "I_xx": 0.001219,
        "mass": 0.175,
    }
    for k, v in d.items():
        np.testing.assert_equal(getattr(e, k), v, k)
    np.testing.assert_allclose(e.diameter, 0.273049106904806)
