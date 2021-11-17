from unittest import TestCase

import numpy as np
import numpy.testing as npt

from frispy import EOM


class TestEOM(TestCase):
    def setUp(self):
        super().setUp()
        # Coordinates for the default
        self.phi = 0
        self.theta = 0
        self.vel = np.array([1, 0, 0])
        self.ang_vel = np.array([0, 0, 62])
        self.kwargs = {
            "area": 0.058556,  # m^2
            "I_zz": 0.002352,  # kg*m^2
            "I_xx": 0.001219,  # kg*m^2
            "mass": 0.175,  # kg
        }

    def test_smoke(self):
        eom = EOM(**self.kwargs)
        assert eom is not None

    def test_eom_has_properties(self):
        eom = EOM(**self.kwargs)
        assert hasattr(eom, "model")
        assert hasattr(eom, "environment")

    def test_compute_forces(self):
        eom = EOM(**self.kwargs)
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        for f in ["F_lift", "F_drag", "F_grav", "F_total", "Acc"]:
            assert f in result
            assert result[f].shape == (3,)
            assert result[f].dtype == np.float

    def test_F_drag_direction(self):
        eom = EOM(**self.kwargs)
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert result["F_drag"][0] < 0  # backwards
        assert result["F_drag"][1] == 0
        assert result["F_drag"][2] == 0

    def test_F_lift_cross_component(self):
        eom = EOM(**self.kwargs)
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert result["F_lift"][1] == 0  # from cross product

    def test_F_grav_direction(self):
        eom = EOM(**self.kwargs)
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert result["F_grav"][0] == 0
        assert result["F_grav"][1] == 0
        assert result["F_grav"][2] < 0  # downwards

    def test_F_total_Acc_relation(self):
        eom = EOM(**self.kwargs)
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        assert all(result["F_total"] == result["Acc"] * eom.mass)

    def test_compute_torques_smoke(self):
        eom = EOM(**self.kwargs)
        result = eom.compute_forces(self.phi, self.theta, self.vel, self.ang_vel)
        result = eom.compute_torques(self.vel, result)
        assert "torque_amplitude" in result
        assert isinstance(result["torque_amplitude"], float)
        for t in ["T_x_lab", "T_y_lab", "T_x", "T_y", "T_z", "T"]:
            assert t in result
            assert result[t].shape == (3,)
            assert result[t].dtype == np.float

    def test_compute_derivatives_smoke(self):
        eom = EOM(**self.kwargs)
        coords = np.array([0, 0, 1, 10, 0, 0, 0, 0, 0, 0, 0, 62])
        der = eom.compute_derivatives(0, coords)
        assert der.shape == (12,)
        assert der.dtype == np.float

    def test_rotation_matrix(self):
        def trig_functions(phi, theta):
            return np.sin(phi), np.cos(phi), np.sin(theta), np.cos(theta)

        r = np.eye(3)  # identity -- no rotation
        assert np.all(EOM.rotation_matrix(*trig_functions(0, 0)) == r)
        # 90 degrees counter clockwise around the primary "x" axis
        r = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        npt.assert_allclose(
            EOM.rotation_matrix(*trig_functions(np.pi / 2, 0)), r, atol=1e-15
        )
        # 90 degrees CCW around the secondary "y" axis
        r = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        npt.assert_allclose(
            EOM.rotation_matrix(*trig_functions(0, np.pi / 2)), r, atol=1e-15
        )
        # 90 degrees CCW around the primary "x" axis then
        # 90 degrees CCW around the secondary "y" axis
        # This permutes the coordinates once
        r = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        npt.assert_allclose(
            EOM.rotation_matrix(*trig_functions(np.pi / 2, np.pi / 2)),
            r,
            atol=1e-15,
        )

    def test_geometric_quantities_case1(self):
        t = EOM(**self.kwargs)
        w = np.array([0, 0, 1])
        res = t.geometric_quantities(0, 0, np.array([1, 0, 0]), w)
        npt.assert_equal(res["w"], w)
        npt.assert_equal(res["w_prime"], w)
        npt.assert_equal(res["w_lab"], w)
        npt.assert_equal(res["rotation_matrix"], np.eye(3))
        assert res["angle_of_attack"] == 0
        npt.assert_equal(res["unit_vectors"]["xhat"], np.array([1, 0, 0]))
        npt.assert_equal(res["unit_vectors"]["yhat"], np.array([0, 1, 0]))
        npt.assert_equal(res["unit_vectors"]["zhat"], np.array([0, 0, 1]))

    def test_geometric_quantities_case2(self):
        t = EOM(**self.kwargs)
        v = np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0])
        w = np.array([0, 0, 1])
        res = t.geometric_quantities(0, 0, v, w)
        npt.assert_equal(res["w"], w)
        npt.assert_equal(res["w_prime"], w)
        npt.assert_equal(res["w_lab"], w)
        npt.assert_equal(res["rotation_matrix"], np.eye(3))
        assert res["angle_of_attack"] == 0
        npt.assert_almost_equal(
            res["unit_vectors"]["xhat"],
            np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0]),
        )
        npt.assert_almost_equal(
            res["unit_vectors"]["yhat"],
            np.array([-1 / np.sqrt(2), 1 / np.sqrt(2), 0]),
        )
        npt.assert_equal(res["unit_vectors"]["zhat"], np.array([0, 0, 1]))

    def test_geometric_quantities_case3(self):
        t = EOM(**self.kwargs)
        theta = np.pi / 4  # 45 degrees
        v = np.array([1, 0, 0])
        w = np.array([0, 0, 1])
        res = t.geometric_quantities(0, theta, v, w)
        npt.assert_almost_equal(res["w"], w)
        npt.assert_equal(res["w_prime"], w)
        npt.assert_almost_equal(
            res["w_lab"], np.array([1 / np.sqrt(2), 0, 1 / np.sqrt(2)])
        )
        npt.assert_almost_equal(
            res["rotation_matrix"],
            np.array(
                [
                    [1 / np.sqrt(2), 0, -1 / np.sqrt(2)],
                    [0, 1, 0],
                    [1 / np.sqrt(2), 0, 1 / np.sqrt(2)],
                ]
            ),
        )
        assert res["angle_of_attack"] == -np.pi / 4
        npt.assert_almost_equal(
            res["unit_vectors"]["xhat"],
            np.array([1 / np.sqrt(2), 0, -1 / np.sqrt(2)]),
        )
        npt.assert_equal(res["unit_vectors"]["yhat"], np.array([0, 1, 0]))
        npt.assert_almost_equal(
            res["unit_vectors"]["zhat"],
            np.array([1 / np.sqrt(2), 0, 1 / np.sqrt(2)]),
        )  # use almost to get around -0 == 0
