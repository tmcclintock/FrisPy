"""
Tests of the ``Disc`` class.
"""

import numpy as np
import numpy.testing as npt

from unittest import TestCase

from frispy import Disc


class TestDisc(TestCase):

    def test_smoke(self):
        d = Disc()
        assert d is not None

    def test_disc_has_trajectory(self):
        d = Disc()
        assert hasattr(d, "_trajectory")
