"""
Tests of the ``Coordinates`` class
"""

from unittest import TestCase

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
