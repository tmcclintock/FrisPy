"""
Tests of the ``Disc`` class.
"""

from unittest import TestCase

from frispy import Disc


class TestDisc(TestCase):
    def test_smoke(self):
        d = Disc()
        assert d is not None

    def test_disc_has_trajectory(self):
        d = Disc()
        assert hasattr(d, "_trajectory")
