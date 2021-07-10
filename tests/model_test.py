"""
Tests of the ``Model`` object.
"""

from unittest import TestCase

import numpy as np
import pytest

from frispy import Model


class TestModel(TestCase):
    def setUp(self):
        super().setUp()
        self.m = Model()

    def test_smoke(self):
        assert self.m is not None

    def test_asserts(self):
        n = "blah"
        v = 0.0
        with pytest.raises(AssertionError):
            Model(blah=v)
        with pytest.raises(AssertionError):
            self.m.set_value(n, v)
        with pytest.raises(AssertionError):
            self.m.set_values({n: v})

    def test_get_value(self):
        n = "PL0"
        v = 1.0
        self.m.set_value(n, v)
        assert self.m.get_value(n) == v

    def test_C_lift(self):
        alphas = np.linspace(-1, 1)
        self.m.set_values({"PL0": 1, "PLa": 1})
        cl = self.m.C_lift(alphas)
        # Linear, strictly increasing
        for i in range(1, len(alphas)):
            assert cl[i] - cl[i - 1] > 0

    def test_C_drag(self):
        alphas = np.linspace(-1, 1, 21)
        self.m.set_values({"PD0": 1, "PDa": 1, "alpha_0": 0})
        cd = self.m.C_drag(alphas)
        # Quadratic, down then up
        for i in range(1, 11):
            assert cd[i] - cd[i - 1] < 0
        for i in range(11, 21):
            assert cd[i] - cd[i - 1] > 0

    def test_C_x(self):
        wx = np.linspace(-1, 1)
        wz = np.linspace(-1, 1)
        self.m.set_values({"PTxwx": 1, "PTxwz": 1})
        cx = self.m.C_x(wx, 0)
        # Linear, stricly increasing
        for i in range(1, len(cx)):
            assert cx[i] - cx[i - 1] > 0
        cx = self.m.C_x(0, wz)
        # Linear, stricly increasing
        for i in range(1, len(cx)):
            assert cx[i] - cx[i - 1] > 0

    def test_C_y(self):
        alphas = np.linspace(-1, 1)
        wy = np.linspace(-1, 1)
        self.m.set_values({"PTy0": 2, "PTywy": 1, "PTya": 1})
        cy = self.m.C_y(alphas, 0)
        # Linear, stricly increasing
        for i in range(1, len(cy)):
            assert cy[i] - cy[i - 1] > 0
        cy = self.m.C_y(0, wy)
        # Linear, stricly increasing
        for i in range(1, len(cy)):
            assert cy[i] - cy[i - 1] > 0
        # Intercept
        assert self.m.C_y(0, 0) == self.m.get_value("PTy0")

    def test_C_z(self):
        wz = np.linspace(-1, 1)
        self.m.set_values({"PTzwz": 1})
        cz = self.m.C_z(wz)
        # Linear, stricly increasing
        for i in range(1, len(cz)):
            assert cz[i] - cz[i - 1] > 0
        # Intercept
        assert self.m.C_z(0) == 0
