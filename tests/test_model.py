"""
Tests of the ``Model`` object.
"""

from unittest import TestCase

import numpy as np

from frispy.model import Model


class TestModel(TestCase):
    def setUp(self):
        super().setUp()
        self.model = Model()

    def test_smoke(self):
        assert self.model is not None

    def test_value_updates(self):
        self.model.PL0 = 1
        assert self.model.PL0 == 1
        self.model.PL0 = 2
        assert self.model.PL0 == 2

    def test_C_lift(self):
        alphas = np.linspace(-1, 1)
        self.model.PL0 = 1
        self.model.PLa = 1
        cl = self.model.C_lift(alphas)
        # Linear, strictly increasing
        for i in range(1, len(alphas)):
            assert cl[i] - cl[i - 1] > 0

    def test_C_drag(self):
        alphas = np.linspace(-1, 1, 21)
        self.model.PD0 = 1
        self.model.PDa = 1
        self.model.alpha_0 = 0
        cd = self.model.C_drag(alphas)
        # Quadratic, down then up
        for i in range(1, 11):
            assert cd[i] - cd[i - 1] < 0
        for i in range(11, 21):
            assert cd[i] - cd[i - 1] > 0

    def test_C_x(self):
        wx = np.linspace(-1, 1)
        wz = np.linspace(-1, 1)
        self.model.PTxwx = 1
        self.model.PTxwz = 1
        cx = self.model.C_x(wx, 0)
        # Linear, stricly increasing
        for i in range(1, len(cx)):
            assert cx[i] - cx[i - 1] > 0
        cx = self.model.C_x(0, wz)
        # Linear, stricly increasing
        for i in range(1, len(cx)):
            assert cx[i] - cx[i - 1] > 0

    def test_C_y(self):
        alphas = np.linspace(-1, 1)
        wy = np.linspace(-1, 1)
        self.model.PTy0 = 1
        self.model.PTywy = 1
        self.model.PTya = 1
        cy = self.model.C_y(alphas, 0)
        # Linear, stricly increasing
        for i in range(1, len(cy)):
            assert cy[i] - cy[i - 1] > 0
        cy = self.model.C_y(0, wy)
        # Linear, stricly increasing
        for i in range(1, len(cy)):
            assert cy[i] - cy[i - 1] > 0
        # Intercept
        assert self.model.C_y(0, 0) == self.model.PTy0

    def test_C_z(self):
        wz = np.linspace(-1, 1)
        self.model.PTzwz = 1
        cz = self.model.C_z(wz)
        # Linear, stricly increasing
        for i in range(1, len(cz)):
            assert cz[i] - cz[i - 1] > 0
        # Intercept
        assert self.model.C_z(0) == 0
