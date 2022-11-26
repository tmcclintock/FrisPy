"""Disc class."""

from typing import Dict, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import OptimizeResult

from frispy.environment import Environment
from frispy.equations_of_motion import EOM
from frispy.model import Model


class Disc:
    """Flying spinning disc object.

    The disc object contains only physical
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

    def __init__(
        self,
        x: float = 0,  # m / s
        y: float = 0,  # m / s
        z: float = 1.0,  # m / s
        vx: float = 10.0,  # m / s
        vy: float = 0,  # m / s
        vz: float = 0,  # m / s
        phi: float = 0,  # rad
        theta: float = 0,  # rad
        gamma: float = 0,  # rad
        dphi: float = 0,  # rad / sec
        dtheta: float = 0,  # rad / sec
        dgamma: float = 62.0,  # rad / sec
        area: float = 0.058556,  # m ^ 2
        I_xx: float = 0.001219,  # kg * m ^ 2
        I_zz: float = 0.002352,  # kg * m ^ 2
        mass: float = 0.175,  # kg
        model: Model = Model(),
        eom: Optional[EOM] = None,
        **kwargs,
    ):
        """Constructor."""
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.phi = phi
        self.theta = theta
        self.gamma = gamma
        self.dphi = dphi
        self.dtheta = dtheta
        self.dgamma = dgamma
        self.area = area
        self.I_xx = I_xx
        self.I_zz = I_zz
        self.mass = mass
        self.model = model
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
        **solver_kwargs,
    ) -> Tuple[Dict[str, np.ndarray], OptimizeResult]:
        """Call the differential equation solver to computethe trajectory.

        The kinematic variables and timesteps are saved
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
            solver_kwargs: extra keyword arguments to pass
                to the :meth:`scipy.integrate.solver_ivp`
        """
        # Pop out these kwargs and take defaults based on our API
        t_span = solver_kwargs.pop("t_span", (0, flight_time))
        t_eval: np.ndarray = solver_kwargs.pop(
            "t_eval", np.linspace(t_span[0], t_span[1], n_times)
        )

        # Call the solver
        result = solve_ivp(
            fun=self.eom.compute_derivatives,
            t_span=t_span,
            y0=[
                self.x,
                self.y,
                self.z,
                self.vx,
                self.vy,
                self.vz,
                self.phi,
                self.theta,
                self.gamma,
                self.dphi,
                self.dtheta,
                self.dgamma,
            ],
            t_eval=t_eval,
            **solver_kwargs,
        )

        return (
            {
                "times": result.t,
                "x": result.y[0],
                "y": result.y[1],
                "z": result.y[2],
                "vx": result.y[3],
                "vy": result.y[4],
                "vz": result.y[5],
                "phi": result.y[6],
                "theta": result.y[7],
                "gamma": result.y[8],
                "dphi": result.y[9],
                "dtheta": result.y[10],
                "dgamma": result.y[11],
            },
            result,
        )

    @property
    def environment(self) -> Environment:
        """Pointer to the envronment.

        TODO: remove
        """
        return self.eom.environment
