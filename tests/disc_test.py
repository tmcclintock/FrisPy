from unittest import TestCase

import pytest

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

    def test_ordered_coordinate_names(self):
        d = Disc()
        assert d.ordered_coordinate_names == [
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
            "phi",
            "theta",
            "gamma",
            "dphi",
            "dtheta",
            "dgamma",
        ]

    def test_disc_has_properties(self):
        d = Disc()
        assert hasattr(d, "trajectory_object")
        assert hasattr(d, "model")
        assert hasattr(d, "environment")
        assert hasattr(d, "eom")

    def test_initial_conditions(self):
        d = Disc()
        assert d.initial_conditions == self.ics
        assert d.current_coordinates == self.ics
        assert d.current_results is None
        ics = self.ics.copy()
        ics["x"] = 1.0
        assert ics != self.ics
        d = Disc(initial_conditions=ics)
        assert d.initial_conditions == ics
        assert d.current_coordinates == ics
        assert d.current_results is None

    def test_reset_initial_conditions(self):
        d = Disc()
        d.current_trajectory = "blah"
        d.current_trajectory_time_points = "lol"
        d.reset_initial_conditions()
        assert d.current_results is None

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

    def test_compute_trajectory_assert_raises_flight_time_and_t_span(self):
        d = Disc()
        with pytest.raises(AssertionError):
            d.compute_trajectory(t_span=(0, 4))

    def test_compute_trajectory_basics(self):
        d = Disc()
        result = d.compute_trajectory()
        for x in d.ordered_coordinate_names:
            assert len(result.times) == len(getattr(result, x))

    def test_compute_trajectory_repeatable(self):
        d = Disc()
        result = d.compute_trajectory()
        for x in d.ordered_coordinate_names:
            assert len(result.times) == len(getattr(result, x))
        result2 = d.compute_trajectory()
        assert all(result.times == result2.times)
        for x in d.ordered_coordinate_names:
            assert len(getattr(result, x)) == len(getattr(result2, x))

    def test_compute_trajectory_return_results(self):
        d = Disc()
        result = d.compute_trajectory()
        result2, scipy_results = d.compute_trajectory(return_scipy_results=True)
        for x in d.ordered_coordinate_names:
            assert len(result.times) == len(getattr(result, x))
        result2 = d.compute_trajectory()
        assert all(result.times == result2.times)
        for x in d.ordered_coordinate_names:
            assert len(getattr(result, x)) == len(getattr(result2, x))
        assert "status" in scipy_results
        assert scipy_results.status >= 0  # -1 is failure

    def test_compute_trajectory_t_span_vs_flight_time(self):
        d = Disc()
        result = d.compute_trajectory(flight_time=3)
        result2 = d.compute_trajectory(t_span=(0, 3), flight_time=None)
        assert all(result.times == result2.times)
        for x in d.ordered_coordinate_names:
            assert len(getattr(result, x)) == len(getattr(result2, x))
