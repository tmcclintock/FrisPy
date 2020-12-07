import pytest

from unittest import TestCase

from frispy import Disc


class TestDisc(TestCase):
    def setUp(self):
        super().setUp()
        self.ics = {
            "x": 0,
            "y": 0,
            "z": 1.0,
            "vx": 10.0,
            "vy": 0,
            "vz": 0,
            "phi": 0,
            "theta": 0,
            "gamma": 0,
            "dphi": 0,
            "dtheta": 0,
            "dgamma": 62.0,
        }

    def test_smoke(self):
        d = Disc()
        assert d is not None

    def test_disc_has_properties(self):
        d = Disc()
        assert hasattr(d, "trajectory")
        assert hasattr(d, "model")
        assert hasattr(d, "environment")
        assert hasattr(d, "eom")

    def test_initial_conditions(self):
        d = Disc()
        assert d.initial_conditions == self.ics
        assert d.current_coordinates == self.ics
        assert d.current_trajectory is None
        ics = self.ics.copy()
        ics["x"] = 1.0
        assert ics != self.ics
        d = Disc(initial_conditions=ics)
        assert d.initial_conditions == ics
        assert d.current_coordinates == ics
        assert d.current_trajectory is None

    def test_reset_initial_conditions(self):
        d = Disc()
        d.current_trajectory = "blah"
        d.reset_initial_conditions()
        assert d.current_trajectory is None

    def test_set_default_initial_conditions(self):
        d = Disc()
        assert d.default_initial_conditions == self.ics
        ics = self.ics.copy()
        ics["x"] = 1.0
        d.set_default_initial_conditions(ics)
        assert ics != self.ics
        assert d.default_initial_conditions == ics
        _ = ics.pop("x")
        with pytest.raises(AssertionError):
            d.set_default_initial_conditions(ics)
