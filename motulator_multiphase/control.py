"""
Control of multi-three-phase machine drives.

The common mode of the three-phase sets is controlled by the vector controllers of
three-phase machines, and the differential modes are current-controlled.

"""

from cmath import exp
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from math import isclose, pi
from typing import Any, Protocol

import numpy as np
from motulator.common.control import (
    PWM,
    ComplexPIController,
    ControlSystem,
    DiscreteComplexController,
    PIController,
    TimeSeries,
)
from motulator.common.utils import abc2complex
from motulator.drive.control.sm import SpeedController

from motulator_multiphase._helpers import empty_array, get_value, wrap
from motulator_multiphase.model import MultiThreePhaseDrive
from motulator_multiphase.utils import default_set_angles, dms_matrix

__all__ = [
    "DifferentialModeController",
    "DifferentialModeSignals",
    "ExternalReferences",
    "Measurements",
    "MultiThreePhaseControlSystem",
]

type Sharing = Sequence[float] | np.ndarray


# %%
class Feedbacks(Protocol):
    """Protocol defining the required fields for the feedback signals."""

    w_M: float  # Mechanical angular speed
    w_c: float  # Angular speed of the coordinate system
    theta_c: float  # Angular position of the coordinate system
    u_s: complex  # Stator voltage (controller coordinates) for the control laws
    u_dc: float  # DC-bus voltage
    i_s: complex  # Current of the parallel-connected machine


class References(Protocol):
    """Protocol defining the required fields for the reference signals."""

    T_s: float  # Sampling period for the next control step
    tau_M: float  # Torque reference
    u_s: complex  # Stator voltage reference
    u_s_lim: complex  # Limited voltage reference (set by the control system)
    w_M: float | None  # Speed reference (set by the control system)
    d_abc: Sequence[float]  # Duty ratios (set by the control system)


class VectorController(Protocol):
    """Protocol defining the interface of the vector controllers of motulator."""

    sensorless: bool
    k_comp: float  # Angle-advance factor

    def get_feedback(
        self,
        u_s_ab: complex,
        i_s_ab: complex,
        w_M_meas: float | None,
        theta_M_meas: float | None,
    ) -> Any:
        """Get feedback signals from measurements."""
        ...

    def compute_output(self, tau_M_ref: float, fbk: Any) -> Any:
        """Compute control output from feedback signals."""
        ...

    def update(self, ref: Any, fbk: Any) -> None:
        """Update controller states."""
        ...

    def post_process(self, ts: TimeSeries) -> None:
        """Post-process controller outputs."""
        ...


@dataclass
class Measurements:
    """Measured signals."""

    i_c_ab: np.ndarray  # Currents of the sets in their own coordinates
    u_dc: float
    w_M: float | None = None
    theta_M: float | None = None


@dataclass
class ExternalReferences:
    """External reference signals."""

    w_M: float | Callable[[float], float] | None = None
    tau_M: float | Callable[[float], float] | None = None
    k_share: Sharing | Callable[[float], Sharing] | None = None


@dataclass
class DifferentialModeSignals:
    """DM signals in the controller coordinates, each of shape (n_sets - 1,)."""

    i: np.ndarray  # Currents
    u: np.ndarray  # Realized voltages
    i_ref: np.ndarray = field(default_factory=empty_array)
    u_ref: np.ndarray = field(default_factory=empty_array)
    u_lim: np.ndarray = field(default_factory=empty_array)  # Limited u_ref


# %%
class DifferentialModeController:
    """
    Current controller of the differential modes.

    Each DM is a resistive-inductive circuit with the leakage inductance `L_ls`. The DM
    currents are controlled in the coordinates of the vector controller using the
    flux-linkage-based complex-vector 2DOF controllers, in the same manner as in the
    current controllers of motulator. The DM current references correspond to the
    current-sharing coefficients of the sets, i.e., the reference of the set k is
    `k_share[k]*i_cm`, where `i_cm` is the common-mode current.

    Parameters
    ----------
    L_ls : float
        Leakage inductance (H) of a set.
    n_sets : int, optional
        Number of three-phase sets, defaults to 2.
    alpha_c : float, optional
        Current-control bandwidth (rad/s), defaults to 2*pi*200.
    alpha_i : float, optional
        Integral-action bandwidth (rad/s), defaults to `alpha_c`. Not used in the
        direct discrete-time design.
    T_s : float, optional
        Sampling period (s). If given, the direct discrete-time design is used instead
        of the continuous-time design. The vector controller must then be a direct
        discrete-time design as well (`k_comp = 0`), since the angle advance of the
        voltage references is common to all modes, and its sampling period must equal
        this value.

    """

    def __init__(
        self,
        L_ls: float,
        n_sets: int = 2,
        alpha_c: float = 2 * pi * 200,
        alpha_i: float | None = None,
        T_s: float | None = None,
    ) -> None:
        self.L_ls = L_ls
        self.n_sets = n_sets
        self.T_s = T_s
        self.discrete = T_s is not None
        self.ctrls: list[ComplexPIController] | list[DiscreteComplexController]
        if T_s is None:
            alpha_i = alpha_c if alpha_i is None else alpha_i
            k_p, k_i = alpha_c + alpha_i, alpha_c * alpha_i
            self.ctrls = [
                ComplexPIController(k_p, k_i, alpha_c) for _ in range(n_sets - 1)
            ]
        else:
            self.ctrls = [
                DiscreteComplexController(alpha_c, T_s) for _ in range(n_sets - 1)
            ]
        self.T = dms_matrix(n_sets)
        self.c_dm = np.zeros(n_sets - 1)  # DM current references per CM current

    def set_sharing(self, k_share: Sharing) -> None:
        """Set the current-sharing coefficients, summing up to the number of sets."""
        k_share = np.asarray(k_share, dtype=float)
        if not np.isclose(np.sum(k_share), self.n_sets):
            raise ValueError(
                "The sharing coefficients must sum up to the number of sets"
            )
        self.c_dm = self.T[1:] @ k_share

    def compute_output(
        self, i_cm: complex, i_dm: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute the DM current references and voltage references."""
        i_dm_ref = self.c_dm * i_cm
        u_dm_ref = [
            ctrl.compute_output(self.L_ls * complex(i_ref), self.L_ls * complex(i))
            for ctrl, i_ref, i in zip(self.ctrls, i_dm_ref, i_dm, strict=True)
        ]
        return i_dm_ref, np.array(u_dm_ref, dtype=complex)

    def update(self, T_s: float, dm: DifferentialModeSignals, w_c: float) -> None:
        """Update the integral states using the realized or limited DM voltages."""
        if self.T_s is not None and not isclose(T_s, self.T_s):
            raise ValueError(
                "The sampling period must equal the design value of the direct "
                "discrete-time DM controller"
            )
        u_dm = dm.u_lim if self.discrete else dm.u
        for ctrl, u in zip(self.ctrls, u_dm, strict=True):
            ctrl.update(T_s, complex(u), w_c)


# %%
class MultiThreePhaseControlSystem(
    ControlSystem[MultiThreePhaseDrive, Measurements, References, Feedbacks]
):
    """
    Vector control system of multi-three-phase machine drives.

    The control is based on the decoupled multi-stator (DMS) transformation, which maps
    the three-phase sets to the common mode (CM) and the differential modes (DMs)
    [#Rub2021]_. The CM is controlled by the given vector controller of a three-phase
    machine, see :func:`~motulator_multiphase.utils.to_parallel`. Its feedback signals
    are the CM voltage and the current of the equivalent machine, i.e., the sum of the
    set currents. The DMs are controlled by the DM current controller in the
    coordinates of the vector controller. Each set has its own PWM. As in
    `VectorControlSystem`, the angles of the voltage references are advanced by
    `k_comp*T_s*w_c`, where `k_comp` is given by the vector controller, and the speed
    controller and the external references work in the same manner. The DM signals are
    logged in the signal group `dm`, see :class:`DifferentialModeSignals`, with the
    shape (N, n_sets - 1). In the post-processed time series, they are in the same
    coordinates as the signals of the vector controller (e.g., in the estimated rotor
    flux coordinates for induction machines).

    Parameters
    ----------
    vector_ctrl : VectorController
        Vector controller of the CM, configured with the parameters of the
        parallel-connected equivalent machine.
    dm_ctrl : DifferentialModeController
        Current controller of the DMs.
    speed_ctrl : SpeedController | PIController | None, optional
        Speed controller. If None (default), torque-control mode is used.
    theta_sets : Sequence[float], optional
        Angles (rad) of the phases a of the sets, defaults to the asymmetrical
        arrangement.
    pwms : Sequence[PWM], optional
        Pulse-width modulators of the sets, defaults to `PWM()` for each set.

    References
    ----------
    .. [#Rub2021] Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation
       for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power
       Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083

    """

    def __init__(
        self,
        vector_ctrl: VectorController,
        dm_ctrl: DifferentialModeController,
        speed_ctrl: SpeedController | PIController | None = None,
        theta_sets: Sequence[float] | np.ndarray | None = None,
        pwms: Sequence[PWM] | None = None,
    ) -> None:
        super().__init__()
        n_sets = dm_ctrl.n_sets
        if theta_sets is None:
            theta_sets = default_set_angles(n_sets)
        if len(theta_sets) != n_sets:
            raise ValueError("The number of the set angles must equal n_sets")
        if pwms is None:
            pwms = [PWM() for _ in range(n_sets)]
        if len(pwms) != n_sets:
            raise ValueError("The number of the PWMs must equal n_sets")
        if dm_ctrl.discrete and vector_ctrl.k_comp != 0:
            raise ValueError(
                "The direct discrete-time DM controller requires a direct "
                "discrete-time vector controller (k_comp = 0)"
            )
        self.vector_ctrl = vector_ctrl
        self.dm_ctrl = dm_ctrl
        self.speed_ctrl = speed_ctrl
        self.n_sets = n_sets
        # Maps the set vectors in their own coordinates to the CM and the DMs
        self.D = dms_matrix(n_sets) * np.exp(1j * np.asarray(theta_sets))
        self.pwms = list(pwms)
        self.k_comp = vector_ctrl.k_comp
        self.ext_ref = ExternalReferences()
        zeros = np.zeros(n_sets - 1, dtype=complex)
        self.dm = DifferentialModeSignals(i=zeros, u=zeros)

    def set_speed_ref(self, ref_fcn: float | Callable[[float], float]) -> None:
        """
        Set the external speed reference for speed-control mode.

        Parameters
        ----------
        ref_fcn : float | Callable[[float], float]
            Speed reference (mechanical rad/s), constant or a function of time.

        """
        if self.speed_ctrl is None:
            raise ValueError("Speed-control mode requires speed_ctrl")
        self.ext_ref.w_M = ref_fcn

    def set_torque_ref(self, ref_fcn: float | Callable[[float], float]) -> None:
        """
        Set the external torque reference for torque-control mode.

        Parameters
        ----------
        ref_fcn : float | Callable[[float], float]
            Torque reference (Nm), constant or a function of time.

        """
        if self.speed_ctrl is not None:
            raise ValueError("The speed controller gives the torque reference")
        self.ext_ref.tau_M = ref_fcn

    def set_current_sharing(
        self, ref_fcn: Sharing | Callable[[float], Sharing]
    ) -> None:
        """
        Set the current-sharing coefficients.

        Parameters
        ----------
        ref_fcn : Sequence[float] | Callable[[float], Sequence[float]]
            Coefficients of the sets, summing up to the number of sets, constant or a
            function of time. Equal sharing, i.e., all coefficients equal to one, is
            used by default.

        """
        self.ext_ref.k_share = ref_fcn

    def get_measurement(self, mdl: MultiThreePhaseDrive) -> Measurements:
        """Get measurements from the sensors."""
        if mdl.machine.n_sets != self.n_sets:
            raise ValueError(
                "The drive model and control system must have equal n_sets"
            )
        i_c_ab = abc2complex(mdl.machine.meas_currents())
        u_dc = mdl.converter.meas_dc_voltage()
        if self.vector_ctrl.sensorless:
            return Measurements(i_c_ab, u_dc)
        w_M = mdl.mechanics.meas_speed()
        theta_M = wrap(mdl.mechanics.meas_position())
        return Measurements(i_c_ab, u_dc, w_M, theta_M)

    def get_feedback(self, meas: Measurements) -> Feedbacks:
        """Get the feedback signals."""
        # The observers use the voltage of the ongoing sampling period, while the
        # control laws use the voltage averaged over the previous and ongoing periods
        u_c_ab = [
            pwm.get_realized_voltage(i, meas.u_dc, average=False)
            for pwm, i in zip(self.pwms, meas.i_c_ab, strict=True)
        ]
        u_zoh_cm = complex(self.D[0] @ u_c_ab)
        u_c_ab = [
            pwm.get_realized_voltage(i, meas.u_dc)
            for pwm, i in zip(self.pwms, meas.i_c_ab, strict=True)
        ]
        u_modes = self.D @ u_c_ab
        i_modes = self.D @ meas.i_c_ab
        # Parallel-connected equivalent machine
        fbk = self.vector_ctrl.get_feedback(
            u_zoh_cm, complex(self.n_sets * i_modes[0]), meas.w_M, meas.theta_M
        )
        exp_mj_theta_c = exp(-1j * fbk.theta_c)
        fbk.u_s = exp_mj_theta_c * complex(u_modes[0])
        fbk.u_dc = meas.u_dc
        self.dm = DifferentialModeSignals(
            i=exp_mj_theta_c * i_modes[1:], u=exp_mj_theta_c * u_modes[1:]
        )
        return fbk

    def compute_output(self, fbk: Feedbacks) -> References:
        """Compute the references and the duty ratios of the sets."""
        ext_ref = self.ext_ref
        if self.speed_ctrl is not None:  # Speed-control mode
            if ext_ref.w_M is None:
                raise ValueError("Speed reference must be set")
            w_M_ref = get_value(ext_ref.w_M, self.t)
            tau_M_ref = self.speed_ctrl.compute_output(w_M_ref, fbk.w_M)
        else:  # Torque-control mode
            if ext_ref.tau_M is None:
                raise ValueError("Torque reference must be set")
            w_M_ref, tau_M_ref = None, get_value(ext_ref.tau_M, self.t)
        ref = self.vector_ctrl.compute_output(tau_M_ref, fbk)
        ref.w_M = w_M_ref

        # DM current control
        dm = self.dm
        if ext_ref.k_share is not None:
            self.dm_ctrl.set_sharing(get_value(ext_ref.k_share, self.t))
        dm.i_ref, dm.u_ref = self.dm_ctrl.compute_output(fbk.i_s / self.n_sets, dm.i)

        # Angle-advanced voltage references of the sets in their own coordinates
        exp_j_theta_c = exp(1j * fbk.theta_c)
        u_modes = exp_j_theta_c * np.array([ref.u_s, *dm.u_ref])
        u_c_ab = self.n_sets * self.D.conj().T @ u_modes
        u_c_ab *= exp(1j * self.k_comp * ref.T_s * fbk.w_c)
        ref.d_abc = [
            d
            for pwm, u in zip(self.pwms, u_c_ab, strict=True)
            for d in pwm(ref.T_s, complex(u), fbk.u_dc, fbk.w_c)
        ]

        # Limited voltages in the controller coordinates
        u_lim_modes = self.D @ [pwm.limited_voltage for pwm in self.pwms]
        u_lim_modes = u_lim_modes / exp_j_theta_c
        ref.u_s_lim, dm.u_lim = complex(u_lim_modes[0]), u_lim_modes[1:]
        return ref

    def compute_standby_output(self, fbk: Feedbacks) -> References:
        """Compute the controller output while the gates are blocked."""
        # The multi-three-phase converter model does not support gate blocking. The
        # message is a separate statement, so that the method is not documented as an
        # abstract method.
        msg = "Multi-three-phase control does not support gate blocking"
        raise NotImplementedError(msg)

    def get_command(self, ref: References) -> tuple[float, Sequence[float]]:
        """Extract the sampling period and the duty ratios from the references."""
        return ref.T_s, ref.d_abc

    def save(self, t: float, **signal_groups: Any) -> None:
        """Extend the base class method by saving the DM signals."""
        super().save(t, dm=self.dm, **signal_groups)

    def update(self, ref: References, fbk: Feedbacks) -> None:
        """Update the controller states."""
        super().update(ref, fbk)
        self.vector_ctrl.update(ref, fbk)
        if self.speed_ctrl is not None:
            self.speed_ctrl.update(ref.T_s, ref.tau_M)
        self.dm_ctrl.update(ref.T_s, self.dm, fbk.w_c)

    def post_process(self) -> TimeSeries:
        """Extend the post-process method."""
        ts = super().post_process()
        theta_c = np.copy(ts.fbk.theta_c)
        self.vector_ctrl.post_process(ts)
        # The vector controller may rotate its signals to other coordinates, e.g., the
        # induction machine controllers to the estimated rotor flux coordinates, and
        # then updates the angle fbk.theta_c. The DM signals are rotated in the same
        # manner.
        T = np.exp(-1j * (ts.fbk.theta_c - theta_c))[:, np.newaxis]
        for name, value in vars(ts.dm).items():
            setattr(ts.dm, name, T * value)
        return ts
