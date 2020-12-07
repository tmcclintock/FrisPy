import pytest
import numpy.testing as npt

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
        assert d.current_trajectory is None
        assert d.current_trajectory_time_points is None
        ics = self.ics.copy()
        ics["x"] = 1.0
        assert ics != self.ics
        d = Disc(initial_conditions=ics)
        assert d.initial_conditions == ics
        assert d.current_coordinates == ics
        assert d.current_trajectory is None
        assert d.current_trajectory_time_points is None

    def test_reset_initial_conditions(self):
        d = Disc()
        d.current_trajectory = "blah"
        d.current_trajectory_time_points = "lol"
        d.reset_initial_conditions()
        assert d.current_trajectory is None
        assert d.current_trajectory_time_points is None

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
        t, traj = d.compute_trajectory()
        assert len(traj) == len(d.initial_conditions_as_ordered_list)
        assert len(t) == len(traj[0])

    def test_compute_trajectory_repeatable(self):
        d = Disc()
        t, traj = d.compute_trajectory()
        assert len(traj) == len(d.initial_conditions_as_ordered_list)
        assert len(t) == len(traj[0])
        t2, traj2 = d.compute_trajectory()
        assert all(t == t2)
        npt.assert_array_equal(traj, traj2)

    def test_compute_trajectory_return_results(self):
        d = Disc()
        t, traj = d.compute_trajectory()
        t2, traj2, results = d.compute_trajectory(return_full_results=True)
        assert all(t == t2)
        npt.assert_array_equal(traj, traj2)
        assert "status" in results
        assert results.status >= 0  # -1 is failure

    def test_compute_trajectory_t_span_vs_flight_time(self):
        d = Disc()
        t, traj = d.compute_trajectory(flight_time=3)
        t2, traj2 = d.compute_trajectory(t_span=(0, 3), flight_time=None)
        assert all(t == t2)
        npt.assert_array_equal(traj, traj2)
