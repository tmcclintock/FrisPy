from collections import OrderedDict, namedtuple
from numbers import Number
from typing import List, Optional, Set

import numpy as np
from scipy.integrate import solve_ivp

from frispy.environment import Environment
from frispy.equations_of_motion import EOM
from frispy.model import Model


class Disc:
    """Flying spinning disc object. The disc object contains only physical
    parameters of the disc and environment that it exists (e.g. gravitational
    acceleration and air density). Note that the default area, mass, and
    inertial moments are for Discraft Ultrastars (175 grams or 0.175 kg).

    All masses are kg, lengths are meters (m) and times are seconds (s). That
    is, these files all use `mks` units. Angular units use radians (rad), and
    angular velocities are in rad/s.

    Args:
        model (Model, optional):
        eom (EOM, optional): the equations of motion
        kwargs: keyword arguments of a numeric type to specify the initial
            conditions of the disc. For example ``x=3`` or ``vz=10.``.
    """

    _default_initial_conditions = OrderedDict(
        {
            "x": 0,
            "y": 0,
            "z": 1.0,
            "vx": 10.0,
            "vy": 0,
            "vz": 0,
            "phi": 0,
            "theta": 0,
            "gamma": 0,
            "dphi": 0,
            "dtheta": 0,
            "dgamma": 62.0,
        }
    )

    _default_physical_attributes = {
        "area": 0.058556,  # m^2
        "I_zz": 0.002352,  # kg*m^2
        "I_xx": 0.001219,  # kg*m^2
        "mass": 0.175,  # kg
    }

    def __init__(
        self, model: Model = Model(), eom: Optional[EOM] = None, **kwargs
    ):
        self.model = model
        self.set_physical_attributes(**kwargs)
        self.set_default_initial_conditions(**kwargs)
        self.reset_initial_conditions()
        self.eom = eom or EOM(
            model=self.model,
            area=self.area,
            I_xx=self.I_xx,
            I_zz=self.I_zz,
            mass=self.mass,
        )

    def compute_trajectory(
        self,
        flight_time: float = 3.0,
        n_times: int = 100,
        return_scipy_results: bool = False,
        **kwargs,
    ):
        """Call the differential equation solver to compute
        the trajectory. The kinematic variables and timesteps are saved
        as the `current_trajectory` attribute, which is a dictionary,
        which is also returned by this function.

        See `the scipy docs
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_
        for more information on the solver.

        .. warning::

           You cannot pass a `flight_time` if `t_span` is a key in
           `solver_args`.

        Args:
            flight_time (float, optional): time in seconds that the simulation
                will run over. Default is 3 seconds.
            n_times (int, optional): default 100. Number of samples in time you
                would like the trajectory. Samples are spaced evenly in time
                from ``(0, flight_time)``.
            return_scipy_results (bool, optional): Default is `False`. Flag to
                indicate whether to return the full results object of the solver.
                See the scipy docs for more information.
            kwargs: extra keyword arguments to pass
                to the :meth:`scipy.integrate.solver_ivp`
        """

        t_span = kwargs.pop("t_span", (0, flight_time))
        t_eval: np.ndarray = kwargs.pop(
            "t_eval", np.linspace(t_span[0], t_span[1], n_times)
        )

        result = solve_ivp(
            fun=self.eom.compute_derivatives,
            t_span=t_span,
            y0=list(self.initial_conditions.values()),
            t_eval=t_eval,
            **kwargs,
        )
        if kwargs.get("dense_output", False):
            return result

        # Set the current coordinates to the last point
        self.current_coordinates = result.y[:, -1]

        # Create the results object
        fpr = Result(times=result.t, *result.y)

        # If specified, return a results object
        if return_scipy_results:
            return fpr, result
        else:
            return fpr

    def reset_initial_conditions(self) -> None:
        """
        Set the initial_conditions of the disc to the default and
        clear the trajectory.
        """
        self.initial_conditions = self.default_initial_conditions
        return

    def set_default_initial_conditions(self, **kwargs) -> None:
        initial_conditions = self._default_initial_conditions.copy()
        valid_keys: Set[str] = set(initial_conditions.keys()).union(
            set(self._default_physical_attributes.keys())
        )
        for key, value in kwargs.items():
            assert key in valid_keys, f"invalid key {key}"
            if key in self._default_physical_attributes:
                pass
            msg = f"invalid type for {key}={value}; {type(value)}"
            assert isinstance(value, Number), msg
            initial_conditions[key] = value
        self.default_initial_conditions = initial_conditions
        return

    def set_physical_attributes(self, **kwargs) -> None:
        for key, value in self._default_physical_attributes.items():
            setattr(self, key, kwargs.get(key, value))
        return

    @property
    def environment(self) -> Environment:
        return self.eom.environment

    @property
    def coordinate_names(self) -> List[str]:
        """
        Names of the kinematic variables
        """
        return list(self._default_initial_conditions.keys())


class Result(
    namedtuple(
        "Result", list(Disc._default_initial_conditions.keys()) + ["times"]
    )
):
    """
    A ``namedtuple`` subclass that contains the coordinate variables
    and a ``times`` attribute. One can reference the variables in the result
    as an attribute ``result.x``.
    """

    pass
