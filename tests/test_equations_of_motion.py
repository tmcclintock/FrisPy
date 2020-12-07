import numpy as np

from unittest import TestCase

from frispy import EOM


class TestEOM(TestCase):
    def setUp(self):
        super().setUp()
        # Coordinates for the default
        self.phi = 0
        self.theta = 0
        self.vel = np.array([1, 0, 0])
        self.ang_vel = np.array([0, 0, 62])

    def test_smoke(self):
        eom = EOM()
        assert eom is not None

    def test_eom_has_properties(self):
        eom = EOM()
        assert hasattr(eom, "trajectory")
        assert hasattr(eom, "model")
        assert hasattr(eom, "environment")

    def test_compute_forces(self):
        eom = EOM()
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        for f in ["F_lift", "F_drag", "F_grav", "F_total", "Acc"]:
            assert f in result
            assert result[f].shape == (3,)
            assert result[f].dtype == np.float

    def test_F_drag_direction(self):
        eom = EOM()
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert result["F_drag"][0] < 0  # backwards
        assert result["F_drag"][1] == 0
        assert result["F_drag"][2] == 0

    def test_F_lift_cross_component(self):
        eom = EOM()
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert result["F_lift"][1] == 0  # from cross product

    def test_F_grav_direction(self):
        eom = EOM()
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert result["F_grav"][0] == 0
        assert result["F_grav"][1] == 0
        assert result["F_grav"][2] < 0  # downwards

    def test_F_total_Acc_relation(self):
        eom = EOM()
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert all(result["F_total"] == result["Acc"] * eom.environment["mass"])

    def test_compute_torques_smoke(self):
        eom = EOM()
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        result = eom.compute_torques(self.vel, self.ang_vel, result)
        assert "torque_amplitude" in result
        assert isinstance(result["torque_amplitude"], float)
        for t in ["T_x_lab", "T_y_lab", "T_x", "T_y", "T_z", "T"]:
            assert t in result
            assert result[t].shape == (3,)
            assert result[t].dtype == np.float

    def test_compute_derivatives_smoke(self):
        eom = EOM()
        coords = [0, 0, 1, 10, 0, 0, 0, 0, 0, 0, 0, 62]
        der = eom.compute_derivatives(coords)
        assert der.shape == (12,)
        assert der.dtype == np.float
