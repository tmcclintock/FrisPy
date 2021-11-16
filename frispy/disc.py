from collections import OrderedDict
from numbers import Number
from typing import Dict, List, Optional

from scipy.integrate import solve_ivp

from frispy.environment import Environment
from frispy.equations_of_motion import EOM
from frispy.model import Model
from frispy.trajectory import Trajectory


class Disc:
    """Flying spinning disc object. The disc object contains only physical
    parameters of the disc and environment that it exists (e.g. gravitational
    acceleration and air density). Note that the default area, mass, and
    inertial moments are for Discraft Ultrastars (175 grams or 0.175 kg).

    All masses are kg, lengths are meters (m) and times are seconds (s). That
    is, these files all use `mks` units. Angular units use radians (rad), and
    angular velocities are in rad/s.

    Args:
      eom (EOM, optional): the equations of motion
      initial_conditions (Dict, optional): initial conditions of the disc in
        flight units are in in "mks" and angles are in radians. By defualt,
        the initial conditions will be that the disc is 1 meter off the
        ground (`z=1`), moving at 10 meters/sec in the `x` direction
        (`vx=10`) and is spinning about it's vertical axis at approximately
        10 revolutions per second (approx. 62 rad/sec, or `dgamma=62`).
        All other kinematic variables are set to 0. This configuration
        results in an angle of attack of 0.
    """

    def __init__(
        self, model: Model = Model(), eom: Optional[EOM] = None, **kwargs
    ):
        self.model = model
        self.eom = eom or EOM(model=self.model)
        self.set_default_initial_conditions(**kwargs)
        self.reset_initial_conditions()

    def compute_trajectory(
        self,
        flight_time: float = 3.0,
        return_scipy_results: bool = False,
        **kwargs,
    ):
        """Call the differential equation solver to compute
        the trajectory. The kinematic variables and timesteps are saved
        as the `current_trajectory` attribute, which is a dictionary,
        which is also returned by this function.

        See `these scipy docs <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_
        for more information on the solver.

        .. todo::

           Type this correctly.

        .. warning::

           You cannot pass a `flight_time` if `t_span` is a key in
           `solver_args`.

        Args:
          flight_time (float, optional): time in seconds that the simulation
            will run over. Default is 3 seconds.
          return_scipy_results (bool, optional): Default is `False`. Flag to
            indicate whether to return the full results object of the solver.
            See the scipy docs for more information.
          solver_args (Dict[str, Any]): extra arguments to pass
            to the :meth:`scipy.integrate.solver_ivp` method used to solve
            the differential equation.
        """
        if "t_span" in kwargs:
            assert (
                flight_time is None
            ), "cannot have t_span in solver_kwargs if flight_time is not None"
            t_span = kwargs.pop("t_span")
        else:
            t_span = (0, flight_time)

        result = solve_ivp(
            fun=self.eom.compute_derivatives,
            t_span=t_span,
            y0=self.initial_conditions_as_ordered_list,
            **kwargs,
        )
        if kwargs.get("dense_output", False):
            return result

        # Set the current coordinates to the last point
        self.current_coordinates = result.y[:, -1]

        # Create the results object
        fpr = {key: result.y[key] for key in self.initial_conditions.keys()}
        fpr["times"] = result.t

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
        initial_conditions = OrderedDict(
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
        assert set(kwargs.keys()).issubset(set(initial_conditions.keys()))
        for key, value in kwargs.items():
            initial_conditions[key] = value
        self.default_initial_conditions = initial_conditions
        return

    @property
    def initial_conditions_names(self) -> List[str]:
        return list(self.initial_conditions.keys())

    @property
    def list_initial_conditions(self) -> List[Number]:
        return list(self.initial_conditions.values())

    @property
    def environment(self) -> Environment:
        return self.eom.environment

    @property
    def trajectory_object(self) -> Trajectory:
        return self.eom.trajectory
