"""
Tests of the ``Coordinates`` class
"""

from unittest import TestCase

import numpy as np
import numpy.testing as npt

from frispy import Coordinates


class TestCoordinates(TestCase):
    def get_coordinates(self):
        coords = {
            "x": 0,
            "y": 0,
            "z": 1,
            "vx": 10,
            "vy": 0,
            "vz": 0,
            "phi": 0,
            "theta": 0,
            "gamma": 0,
            "dphi": 0,
            "dtheta": 0,
            "dgamma": 50,
        }
        return Coordinates(**coords)

    def test_smoke(self):
        co = self.get_coordinates()
        assert co is not None

    def test_position(self):
        co = self.get_coordinates()
        assert all(co.position == np.array([co.x, co.y, co.z]))

    def test_velocity(self):
        co = self.get_coordinates()
        assert all(co.velocity == np.array([co.vx, co.vy, co.vz]))

    def test_rotation_matrix_0_0(self):
        co = self.get_coordinates()
        co.phi = 0
        co.theta = 0
        # identity -- no rotation
        r = np.eye(3)
        assert np.all(co.rotation_matrix == r)

    def test_rotation_matrix_pi2_0(self):
        co = self.get_coordinates()
        # assert np.all(Trajectory.rotation_matrix(0, 0) == r)
        # 90 degrees counter clockwise around the primary "x" axis
        co.phi = np.pi / 2
        r = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        npt.assert_allclose(co.rotation_matrix, r, atol=1e-15)

    def test_rotation_matrix_0_pi2(self):
        co = self.get_coordinates()
        # 90 degrees CCW around the secondary "y" axis
        r = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        co.phi = 0
        co.theta = np.pi / 2
        npt.assert_allclose(co.rotation_matrix, r, atol=1e-15)

    def test_rotation_matrix_pi2_pi2(self):
        co = self.get_coordinates()
        # 90 degrees CCW around the primary "x" axis then
        # 90 degrees CCW around the secondary "y" axis
        # This permutes the coordinates once
        r = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        co.phi = co.theta = np.pi / 2
        npt.assert_allclose(co.rotation_matrix, r, atol=1e-15)
