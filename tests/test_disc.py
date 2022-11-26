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


def test_initial_conditions():
    d = Disc()
    for key, value in d._default_initial_conditions.items():
        assert d.default_initial_conditions[key] == value


def test_initial_conditions_kwarg():
    d = Disc(vz=2.0)
    for key, value in d._default_initial_conditions.items():
        if key == "vz":
            assert d.default_initial_conditions["vz"] == 2.0
        else:
            assert d.default_initial_conditions[key] == value


def test_physical_attribute_kwarg():
    d = Disc(mass=12345, area=0.1234)
    assert d.mass == 12345
    assert d.area == 0.1234
    assert d.eom.diameter == 2 * np.sqrt(d.area / np.pi)


def test_compute_trajectory_basics():
    d = Disc()
    result = d.compute_trajectory()
    for x in d.coordinate_names:
        assert len(result.times) == len(getattr(result, x))


def test_compute_trajectory_repeatable():
    d = Disc()
    result = d.compute_trajectory()
    for x in d.coordinate_names:
        assert len(result.times) == len(getattr(result, x))
    result2 = d.compute_trajectory()
    assert all(result.times == result2.times)
    for x in d.coordinate_names:
        assert len(getattr(result, x)) == len(getattr(result2, x))


def test_compute_trajectory_return_results():
    d = Disc()
    result = d.compute_trajectory()
    result2, scipy_results = d.compute_trajectory(return_scipy_results=True)
    for x in d.coordinate_names:
        assert len(result.times) == len(getattr(result, x))
    result2 = d.compute_trajectory()
    assert all(result.times == result2.times)
    for x in d.coordinate_names:
        assert len(getattr(result, x)) == len(getattr(result2, x))
    assert "status" in scipy_results
    assert scipy_results.status >= 0  # -1 is failure


def test_compute_trajectory_t_span_vs_flight_time():
    d = Disc()
    result = d.compute_trajectory(flight_time=3)
    result2 = d.compute_trajectory(t_span=(0, 3), flight_time=None)
    assert all(result.times == result2.times)
    for x in d.coordinate_names:
        assert len(getattr(result, x)) == len(getattr(result2, x))
