import FrisPy
import numpy as np
import numpy.testing as npt
import pytest

def test_Trajectory():
    #Smoke test
    t = FrisPy.Trajectory()

def test_trajectory_defaults():
    t = FrisPy.Trajectory()
    allowed_keys = ['x', 'y', 'z', 'vx', 'vy', 'vz',
                    'phi', 'theta', 'gamma',
                    'phidot', 'thetadot', 'gammadot']
    default_values = [0, 0, 1, 10, 0, 0,
                      0, 0, 0, 0, 0, 50]
    for i, (k, v) in enumerate(zip(allowed_keys, default_values)):
        npt.assert_equal(getattr(t, k), default_values[i])

    d = FrisPy.Disc()
    
