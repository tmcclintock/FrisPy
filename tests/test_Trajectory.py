"""
Tests of the ``Trajectory`` class.
"""

import numpy as np
import numpy.testing as npt

from unittest import TestCase

from frispy import Trajectory


class TestTrajectory(TestCase):

    def test_smoke(self):
        t = Trajectory()
        assert t is not None

    def test_trajectory_default_start_values(self):
        t = Trajectory()
        truth = {'x': 0,
                 'y': 0,
                 'z': 1,
                 'vx': 10,
                 'vy':0,
                 'vz':0,
                 "phi":0,
                 "theta":0,
                 "gamma":0,
                 "phidot":0,
                 "thetadot":0,
                 "gammadot":50,
        }
        for k, v in truth.items():
            assert t.initial_conditions[k] == v
        truth_arr = np.array([truth[k] for k in t._coord_order])
        np.testing.assert_equal(t.initial_conditions_array, truth_arr)
    
    def test_rotation_matrix(self):
        r = np.eye(3)  # identity -- no rotation
        assert np.all(Trajectory.rotation_matrix(0, 0) == r)
        # 90 degrees counter clockwise around the primary "x" axis
        r = np.array([
            [1, 0, 0],
            [0, 0, 1],
            [0, -1, 0]
        ])
        np.testing.assert_allclose(
            Trajectory.rotation_matrix(np.pi / 2, 0),
            r,
            atol=1e-15
        )
        # 90 degrees CCW around the secondary "y" axis
        r = np.array([
            [0, 0, -1],
            [0, 1, 0],
            [1, 0, 0]
        ])
        np.testing.assert_allclose(
            Trajectory.rotation_matrix(0, np.pi / 2),
            r,
            atol=1e-15
        )
        # 90 degrees CCW around the primary "x" axis then
        # 90 degrees CCW around the secondary "y" axis
        # This permutes the coordinates once
        r = np.array([
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 0]
        ])
        np.testing.assert_allclose(
            Trajectory.rotation_matrix(np.pi / 2, np.pi / 2),
            r,
            atol=1e-15
        )
