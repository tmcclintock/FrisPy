import numpy as np

from frispy.disc import Disc


def test_smoke():
    d = Disc()
    assert d is not None


def test_disc_has_properties():
    d = Disc()
    assert hasattr(d, "model")
    assert hasattr(d, "environment")
    assert hasattr(d, "eom")


def test_physical_attribute_kwarg():
    d = Disc(mass=12345, area=0.1234)
    assert d.mass == 12345
    assert d.area == 0.1234
    assert d.eom.diameter == 2 * np.sqrt(d.area / np.pi)


def test_compute_trajectory_basics():
    d = Disc()
    dct, _ = d.compute_trajectory()
    assert len(dct) == 13
    n_times = len(dct["times"])
    for v in dct.values():
        assert len(v) == n_times


def test_compute_trajectory_repeatable():
    d = Disc()
    dct1, _ = d.compute_trajectory()
    dct2, _ = d.compute_trajectory()
    for k, v in dct1.items():
        np.testing.assert_array_equal(dct2[k], v)


def test_compute_trajectory_t_span_vs_flight_time():
    d = Disc()
    dct1, _ = d.compute_trajectory(flight_time=3)
    dct2, _ = d.compute_trajectory(t_span=(0, 3), flight_time=None)
    np.testing.assert_array_equal(dct1["times"], dct2["times"])
    for k, v in dct1.items():
        np.testing.assert_array_equal(dct2[k], v)
