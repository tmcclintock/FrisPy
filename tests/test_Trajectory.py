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
            assert getattr(t, k) == v
    
