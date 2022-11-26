"""Disc class."""

from dataclasses import dataclass
from typing import Dict, Tuple, Type

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import OptimizeResult

from frispy.equations_of_motion import EOM
from frispy.model import Model


@dataclass
class Disc:
    """Flying spinning disc object.

    One can specify initial coordinates, initial kinematic variables, the two
    non-zero elements of the moments of inertia tensor, the face area, mass,
    air density, and gravitational acceleration. Note that the default area,
    mass, and inertial moments are for Discraft Ultrastars (0.175 kg).

    MKS units used throughout. Angular units use radians (rad), and
    angular velocities are in rad/s.

    Args:
        x: horizontal position
        y: horizontal position
        z: vertical position
        vx: horizontal speed
        vy: horizontal speed
        vz: vertical speed
        phi: first Euler angle. Pitch angle when velocity is horizontal
        theta: second Euler angle. Roll angle when velocity is horizontal
        gamma: third Euler angle. Spin angle
        dphi: angular velocity
        dtheta: angular velocity
        dgamma: spin angular velocity
        area: face area
        I_xx: pitch and roll moments of inertia
        I_zz: spin moment of inertia
        air_density: default is sea level
        g: gravitational acceleration
        model: defines how accelerations depend on velocities
        eom: equations of motion
    """

    x: float = 0  # m / s
    y: float = 0  # m / s
    z: float = 1.0  # m / s
    vx: float = 10.0  # m / s
    vy: float = 0  # m / s
    vz: float = 0  # m / s
    phi: float = 0  # rad
    theta: float = 0  # rad
    gamma: float = 0  # rad
    dphi: float = 0  # rad / sec
    dtheta: float = 0  # rad / sec
    dgamma: float = 62.0  # rad / sec
    area: float = 0.058556  # m ^ 2
    I_xx: float = 0.001219  # kg * m ^ 2
    I_zz: float = 0.002352  # kg * m ^ 2
    mass: float = 0.175  # kg
    air_density: float = 1.225  # kg / m ^ 3
    g: float = 9.81  # m / s ^ 2
    model: Model = Model()
    eom_class: Type = EOM

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
            flight_time: time in seconds that the simulation
                will run over
            n_times: Number of samples in time for the trajectory.
                Samples are spaced evenly in time from ``(0, flight_time)``.
            solver_kwargs: extra keyword arguments to pass
                to the :meth:`scipy.integrate.solver_ivp`

        Returns:
            first element is a dict of all kinematic variables and the
            second is the ``solver_ivp`` results object
        """
        # Pop out these kwargs and take defaults based on our API
        t_span = solver_kwargs.pop("t_span", (0, flight_time))
        t_eval: np.ndarray = solver_kwargs.pop(
            "t_eval", np.linspace(t_span[0], t_span[1], n_times)
        )

        # Instantiate the equations of motion
        eom = self.eom_class(
            model=self.model,
            area=self.area,
            I_xx=self.I_xx,
            I_zz=self.I_zz,
            mass=self.mass,
            air_density=self.air_density,
            g=self.g,
        )

        # Call the solver
        result = solve_ivp(
            fun=eom.compute_derivatives,
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
