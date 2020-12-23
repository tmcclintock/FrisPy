"""
Tests of the ``Trajectory`` class.
"""

from unittest import TestCase

import numpy as np
import numpy.testing as npt

from frispy import Trajectory


class TestTrajectory(TestCase):
    def test_smoke(self):
        t = Trajectory()
        assert t is not None

    def test_trajectory_default_start_values(self):
        t = Trajectory()
        truth = {
            "x": 0,
            "y": 0,
            "z": 1,
            "vx": 10,
            "vy": 0,
            "vz": 0,
            "phi": 0,
            "theta": 0,
            "gamma": 0,
            "phidot": 0,
            "thetadot": 0,
            "gammadot": 50,
        }
        for k, v in truth.items():
            assert t.initial_conditions[k] == v
        truth_arr = np.array([truth[k] for k in t._coord_order])
        npt.assert_equal(t.initial_conditions_array, truth_arr)

    def test_rotation_matrix(self):
        r = np.eye(3)  # identity -- no rotation
        assert np.all(Trajectory.rotation_matrix(0, 0) == r)
        # 90 degrees counter clockwise around the primary "x" axis
        r = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        npt.assert_allclose(
            Trajectory.rotation_matrix(np.pi / 2, 0), r, atol=1e-15
        )
        # 90 degrees CCW around the secondary "y" axis
        r = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        npt.assert_allclose(
            Trajectory.rotation_matrix(0, np.pi / 2), r, atol=1e-15
        )
        # 90 degrees CCW around the primary "x" axis then
        # 90 degrees CCW around the secondary "y" axis
        # This permutes the coordinates once
        r = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        npt.assert_allclose(
            Trajectory.rotation_matrix(np.pi / 2, np.pi / 2), r, atol=1e-15
        )

    def test_velocity(self):
        t = Trajectory()
        npt.assert_equal(t.velocity, np.array([10.0, 0, 0]))

    def test_calculate_intermediate_quantities_case1(self):
        t = Trajectory()
        phi = 0
        theta = 0
        v = np.array([1, 0, 0])
        w = np.array([0, 0, 1])
        res = t.calculate_intermediate_quantities(phi, theta, v, w)
        npt.assert_equal(res["w"], w)
        npt.assert_equal(res["w_prime"], w)
        npt.assert_equal(res["w_lab"], w)
        npt.assert_equal(res["rotation_matrix"], np.eye(3))
        assert res["angle_of_attack"] == 0
        npt.assert_equal(res["unit_vectors"]["xhat"], np.array([1, 0, 0]))
        npt.assert_equal(res["unit_vectors"]["yhat"], np.array([0, 1, 0]))
        npt.assert_equal(res["unit_vectors"]["zhat"], np.array([0, 0, 1]))

    def test_calculate_intermediate_quantities_case2(self):
        t = Trajectory()
        phi = 0
        theta = 0
        v = np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0])
        w = np.array([0, 0, 1])
        res = t.calculate_intermediate_quantities(phi, theta, v, w)
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

    def test_calculate_intermediate_quantities_case3(self):
        t = Trajectory()
        phi = 0
        theta = np.pi / 4  # 45 degrees
        v = np.array([1, 0, 0])
        w = np.array([0, 0, 1])
        res = t.calculate_intermediate_quantities(phi, theta, v, w)
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
