"""
Continuous-time models for multi-three-phase machine drives.

Peak-valued complex space vectors are used. The voltages and currents of the
three-phase sets are space vectors in the sets' own coordinates. The mechanical
system models, the machine parameter classes, and the simulation class of motulator
are re-exported for convenience.

"""

from collections.abc import Callable, Sequence
from dataclasses import InitVar, dataclass, field
from types import SimpleNamespace
from typing import Any

import numpy as np
from motulator.common.model import Model, Subsystem, SubsystemTimeSeries
from motulator.common.utils import abc2complex, complex2abc
from motulator.drive.model import (
    ExternalRotorSpeed,
    InductionMachine,
    InductionMachineInvGammaPars,
    InductionMachinePars,
    MechanicalSystem,
    SaturatedSynchronousMachinePars,
    Simulation,
    SpatialSaturatedSynchronousMachinePars,
    SynchronousMachine,
    SynchronousMachinePars,
    TwoMassMechanicalSystem,
)

from motulator_multiphase._helpers import empty_array
from motulator_multiphase.utils import (
    MachinePars,
    default_set_angles,
    dms_matrix,
    to_parallel,
)

__all__ = [
    "ExternalRotorSpeed",
    "InductionMachineInvGammaPars",
    "InductionMachinePars",
    "MechanicalSystem",
    "MultiThreePhaseConverter",
    "MultiThreePhaseDrive",
    "MultiThreePhaseMachine",
    "SaturatedSynchronousMachinePars",
    "Simulation",
    "SpatialSaturatedSynchronousMachinePars",
    "SynchronousMachinePars",
    "TwoMassMechanicalSystem",
]


# %%
@dataclass
class ConverterInputs:
    """Input variables."""

    q_c_ab: np.ndarray  # Switching states of the sets, shape (n_sets,)
    s_c_ab: np.ndarray  # Current directions at the latest switching instant
    i_c_ab: np.ndarray  # Currents of the sets, shape (n_sets,)
    enabled: bool = True  # Gate enable, always True (blocking is not supported)


@dataclass
class ConverterOutputs:
    """Output variables for interconnection."""

    u_c_ab: np.ndarray  # Voltages of the sets, shape (n_sets,)
    u_dc: float


class MultiThreePhaseConverter(Subsystem):
    """
    Lossless multi-three-phase converter with constant DC-bus voltage.

    The converter consists of `n_sets` three-phase converters sharing the DC bus. The
    switching states, voltages, and currents of each set are space vectors in the
    set's own coordinates. The leg states are ordered set by set, `[a1, b1, c1, a2,
    ...]`. As in the three-phase converter, the switches are ideal except for the dead
    time: when a leg is blanked, its state is `(1 - sign(i))/2`. Gate blocking is not
    supported, see `set_enable`.

    Parameters
    ----------
    u_dc : float
        DC-bus voltage (V).
    n_sets : int, optional
        Number of three-phase sets, defaults to 2.
    t_d : float, optional
        Dead time (s), defaults to 0.
    sign : Callable[[np.ndarray], np.ndarray], optional
        Function of the phase currents (A) determining the leg states during blanking,
        defaults to `np.sign`.

    """

    def __init__(
        self,
        u_dc: float,
        n_sets: int = 2,
        t_d: float = 0.0,
        sign: Callable[[np.ndarray], np.ndarray] = np.sign,
    ) -> None:
        self.u_dc = u_dc
        self.n_sets = n_sets
        self.n_legs = 3 * n_sets
        self.t_d = t_d
        self.sign = sign
        self.inp: ConverterInputs = ConverterInputs(
            q_c_ab=np.zeros(n_sets, dtype=complex),
            s_c_ab=np.zeros(n_sets, dtype=complex),
            i_c_ab=np.zeros(n_sets, dtype=complex),
        )
        self.out: ConverterOutputs = ConverterOutputs(
            u_c_ab=np.zeros(n_sets, dtype=complex), u_dc=u_dc
        )
        self.state = None
        self.ac_load: Subsystem | None = None  # No AC port without gate blocking
        self.margin = np.inf  # No diode-state transitions without gate blocking

    def set_enable(self, enabled: bool) -> None:
        """Set the gate enable, raising an error if the gates are to be blocked."""
        if not enabled:
            raise NotImplementedError(
                "The multi-three-phase converter does not support gate blocking"
            )

    def update_diodes(self) -> None:
        """Not used, since the gates are never blocked."""
        raise NotImplementedError(
            "The multi-three-phase converter does not support gate blocking"
        )

    def set_gate_signals(self, q_abc: np.ndarray, b_abc: np.ndarray) -> None:
        """Set the switching states based on the gate and blanking signals."""
        s_abc = self.sign(complex2abc(self.inp.i_c_ab).T.ravel())
        q_abc = q_abc + 0.5 * b_abc * (1 - s_abc)
        self.inp.q_c_ab = abc2complex(np.reshape(q_abc, (self.n_sets, 3)).T)
        self.inp.s_c_ab = abc2complex(np.reshape(s_abc, (self.n_sets, 3)).T)

    def set_outputs(self, t: float) -> None:
        """Set output variables."""
        self.out.u_c_ab = self.inp.q_c_ab * self.out.u_dc

    def meas_dc_voltage(self) -> float:
        """Measure the DC-bus voltage (V)."""
        return self.out.u_dc

    def rhs(self, t: float) -> list[complex]:
        """No states."""
        return []

    def create_time_series(
        self, t: np.ndarray
    ) -> tuple[str, "MultiThreePhaseConverterTimeSeries"]:
        """Create time series."""
        return "converter", MultiThreePhaseConverterTimeSeries(t, self)


@dataclass
class MultiThreePhaseConverterTimeSeries(SubsystemTimeSeries):
    """Continuous time series, the set quantities of shape (n_sets, N)."""

    t: InitVar[np.ndarray]
    subsystem: InitVar[MultiThreePhaseConverter]
    u_dc: np.ndarray = field(default_factory=empty_array)
    q_c_ab: np.ndarray = field(default_factory=empty_array)
    s_c_ab: np.ndarray = field(default_factory=empty_array)
    enabled: np.ndarray = field(default_factory=empty_array)
    i_c_ab: np.ndarray = field(default_factory=empty_array)
    u_c_ab: np.ndarray = field(default_factory=empty_array)

    def __post_init__(self, t: np.ndarray, subsystem: MultiThreePhaseConverter) -> None:
        self.u_dc = np.full(np.size(t), subsystem.u_dc)

    def compute_zoh_input_derived_signals(
        self, t: np.ndarray, subsystem: MultiThreePhaseConverter
    ) -> None:
        """Compute the set voltages."""
        # The history has the shape (N, n_sets), or (0,) if it is empty
        self.q_c_ab = np.reshape(self.q_c_ab, (-1, subsystem.n_sets)).T
        self.s_c_ab = np.reshape(self.s_c_ab, (-1, subsystem.n_sets)).T
        self.u_c_ab = self.q_c_ab * self.u_dc


# %%
@dataclass
class MachineInputs:
    """Input variables."""

    u_c_ab: np.ndarray  # Voltages of the sets, shape (n_sets,)
    w_M: float = 0.0


@dataclass
class MachineOutputs:
    """Output variables for interconnection."""

    i_c_ab: np.ndarray  # Currents of the sets, shape (n_sets,)
    tau_M: float


class MultiThreePhaseMachine(Subsystem):
    """
    Multi-three-phase machine model.

    The machine has `n_sets` three-phase winding sets. The model is based on the
    decoupled multi-stator (DMS) transformation, which maps the sets to the common mode
    (CM) and the differential modes (DMs) [#Rub2021]_, see :func:`dms_matrix`. The CM is
    equivalent to a three-phase machine whose sets are connected in parallel. It is
    modeled by the three-phase machine model `cm` (`SynchronousMachine` or
    `InductionMachine`) with the parameters given by :func:`to_parallel`. Each DM is a
    resistive-inductive circuit, consisting of the stator resistance and the leakage
    inductance `L_ls` of a set, which is not coupled to the other sets. The DMs produce
    no torque, and they are modeled in stationary coordinates. The states are those of
    the CM machine, followed by the DM flux linkages.

    Parameters
    ----------
    par : MachinePars
        Parameters of one three-phase set, when all sets carry equal currents.
    L_ls : float
        Leakage inductance (H) of a set, which is not coupled to the other sets.
    n_sets : int, optional
        Number of three-phase sets, defaults to 2.
    theta_sets : Sequence[float], optional
        Angles (rad) of the phases a of the sets, defaults to the asymmetrical
        arrangement, see :func:`default_set_angles`.

    Notes
    -----
    The resistance of the DMs is the stator resistance `par.R_s` of a set, which is the
    same in the Γ and inverse-Γ models of the induction machine.

    References
    ----------
    .. [#Rub2021] Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation
       for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power
       Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083

    """

    def __init__(
        self,
        par: MachinePars,
        L_ls: float,
        n_sets: int = 2,
        theta_sets: Sequence[float] | np.ndarray | None = None,
    ) -> None:
        if theta_sets is None:
            theta_sets = default_set_angles(n_sets)
        if len(theta_sets) != n_sets:
            raise ValueError("The number of the set angles must equal n_sets")
        self.n_sets = n_sets
        self.theta_sets = np.asarray(theta_sets, dtype=float)
        self.R_s = par.R_s
        self.L_ls = L_ls
        # Maps the set vectors in their own coordinates to the CM and the DMs
        self.D = dms_matrix(n_sets) * np.exp(1j * self.theta_sets)
        par_cm = to_parallel(par, n_sets)
        self.cm: SynchronousMachine | InductionMachine
        if isinstance(par_cm, (InductionMachinePars, InductionMachineInvGammaPars)):
            self.cm = InductionMachine(par_cm)
        else:
            self.cm = SynchronousMachine(par_cm)
        self._dm_fields = [f"psi_dm{u}_ab" for u in range(1, n_sets)]
        self.state = SimpleNamespace(
            **vars(self.cm.state), **dict.fromkeys(self._dm_fields, 0j)
        )
        self.inp: MachineInputs = MachineInputs(u_c_ab=np.zeros(n_sets, dtype=complex))
        self.out: MachineOutputs = MachineOutputs(
            i_c_ab=np.zeros(n_sets, dtype=complex), tau_M=0.0
        )
        self.set_outputs(0.0)

    @property
    def psi_dm_ab(self) -> np.ndarray:
        """DM flux linkages (Vs), shape (n_sets - 1,)."""
        return np.array(
            [getattr(self.state, f) for f in self._dm_fields], dtype=complex
        )

    def set_states(self, state_list: list[complex], index: int) -> int:
        """Set the states, also those of the CM machine."""
        self.cm.set_states(state_list, index)
        index = super().set_states(state_list, index)
        # The CM machine wraps its rotor angle, and the next integration interval
        # starts from the wrapped angle
        for name, value in vars(self.cm.state).items():
            setattr(self.state, name, value)
        return index

    def init_state_history(self) -> None:
        """Initialize the state history, the CM states in the CM machine."""
        self.cm.init_state_history()
        self._history = {f: [] for f in self._dm_fields}

    def extend_state_history(self, sol_y: list[list[float]], index: int) -> int:
        """Extend the state history, the CM states in the CM machine."""
        index = self.cm.extend_state_history(sol_y, index)
        return super().extend_state_history(sol_y, index)

    def compute_currents(self, i_s_ab: Any, psi_dm_ab: Any) -> np.ndarray:
        """
        Compute the currents of the sets.

        Parameters
        ----------
        i_s_ab : complex | ndarray, shape (N,)
            Current (A) of the parallel-connected machine, i.e., `n_sets` times the CM.
        psi_dm_ab : ndarray, shape (n_sets - 1,) or (n_sets - 1, N)
            DM flux linkages (Vs).

        Returns
        -------
        ndarray, shape (n_sets,) or (n_sets, N)
            Currents (A) of the sets in their own coordinates.

        """
        i_modes = np.array([i_s_ab / self.n_sets, *(psi_dm_ab / self.L_ls)])
        return self.n_sets * self.D.conj().T @ i_modes

    def _set_cm_inputs(self) -> None:
        """Pass the CM voltage and the speed to the CM machine."""
        self.cm.inp.u_s_ab = complex(self.D[0] @ self.inp.u_c_ab)
        self.cm.inp.w_M = self.inp.w_M

    def set_outputs(self, t: float) -> None:
        """Set output variables."""
        self._set_cm_inputs()
        self.cm.set_outputs(t)
        self.out.i_c_ab = self.compute_currents(self.cm.out.i_s_ab, self.psi_dm_ab)
        self.out.tau_M = self.cm.out.tau_M

    def rhs(self, t: float) -> list[complex]:
        """Compute state derivatives."""
        # The mechanics may have updated the speed after set_outputs
        self._set_cm_inputs()
        u_dm_ab = self.D[1:] @ self.inp.u_c_ab
        d_psi_dm_ab = u_dm_ab - self.R_s * self.psi_dm_ab / self.L_ls
        return [*self.cm.rhs(t), *d_psi_dm_ab]

    def meas_currents(self) -> np.ndarray:
        """Measure the phase currents (A) of the sets, shape (3, n_sets)."""
        return complex2abc(self.out.i_c_ab)

    def create_time_series(
        self, t: np.ndarray
    ) -> tuple[str, "MultiThreePhaseMachineTimeSeries"]:
        """Create time series."""
        return "machine", MultiThreePhaseMachineTimeSeries(t, self)


@dataclass
class MultiThreePhaseMachineTimeSeries(SubsystemTimeSeries):
    """
    Continuous time series of the multi-three-phase machine.

    The time series of the CM machine is `cm`. The stator flux linkage, current, and
    voltage of the parallel-connected machine are also given as `psi_s_ab`, `i_s_ab`,
    and `u_s_ab`, so that the plotting functions of three-phase drives can be used. The
    set quantities have the shape (n_sets, N) and the DM quantities (n_sets - 1, N).

    """

    t: InitVar[np.ndarray]
    subsystem: InitVar[MultiThreePhaseMachine]
    cm: Any = field(init=False)  # Time series of the CM machine
    # States
    psi_dm_ab: np.ndarray = field(default_factory=empty_array)
    # Outputs
    i_c_ab: np.ndarray = field(default_factory=empty_array)
    tau_M: np.ndarray = field(default_factory=empty_array)
    # Inputs
    u_c_ab: np.ndarray = field(default_factory=empty_array)
    w_M: np.ndarray = field(default_factory=empty_array)
    # Derived signals
    i_dm_ab: np.ndarray = field(default_factory=empty_array)
    psi_s_ab: np.ndarray = field(default_factory=empty_array)
    i_s_ab: np.ndarray = field(default_factory=empty_array)
    u_s_ab: np.ndarray = field(default_factory=empty_array)

    def __post_init__(self, t: np.ndarray, subsystem: MultiThreePhaseMachine) -> None:
        _, self.cm = subsystem.cm.create_time_series(t)
        self.psi_dm_ab = np.reshape(
            [subsystem._history[f] for f in subsystem._dm_fields],
            (subsystem.n_sets - 1, np.size(t)),
        )
        self.i_dm_ab = self.psi_dm_ab / subsystem.L_ls
        self.psi_s_ab, self.tau_M = self.cm.psi_s_ab, self.cm.tau_M
        self._set_currents(subsystem)

    def _set_currents(self, subsystem: MultiThreePhaseMachine) -> None:
        self.i_s_ab = self.cm.i_s_ab
        self.i_c_ab = subsystem.compute_currents(self.i_s_ab, self.psi_dm_ab)

    def compute_input_derived_signals(
        self, t: np.ndarray, subsystem: MultiThreePhaseMachine
    ) -> None:
        """Pass the inputs to the CM machine and update the currents."""
        self.u_s_ab = subsystem.D[0] @ self.u_c_ab
        self.cm.u_s_ab, self.cm.w_M = self.u_s_ab, self.w_M
        self.cm.compute_input_derived_signals(t, subsystem.cm)
        self._set_currents(subsystem)


# %%
class MultiThreePhaseDrive(Model[MultiThreePhaseConverter]):
    """
    Continuous-time model of a multi-three-phase machine drive.

    Parameters
    ----------
    machine : MultiThreePhaseMachine
        Machine model.
    mechanics : MechanicalSystem | TwoMassMechanicalSystem | ExternalRotorSpeed
        Mechanical system model.
    converter : MultiThreePhaseConverter
        Converter model.
    pwm : bool, optional
        Enable the carrier comparison, defaults to False.
    delay : int, optional
        Computational delay (samples), defaults to 1.

    """

    def __init__(
        self,
        machine: MultiThreePhaseMachine,
        mechanics: MechanicalSystem | TwoMassMechanicalSystem | ExternalRotorSpeed,
        converter: MultiThreePhaseConverter,
        pwm: bool = False,
        delay: int = 1,
    ) -> None:
        if machine.n_sets != converter.n_sets:
            raise ValueError("The machine and converter must have equal n_sets")
        self.machine = machine
        self.mechanics = mechanics
        subsystems = [converter, machine, mechanics]
        connections = {
            (machine, "u_c_ab"): (converter, "u_c_ab"),
            (machine, "w_M"): (mechanics, "w_M"),
            (converter, "i_c_ab"): (machine, "i_c_ab"),
            (mechanics, "tau_M"): (machine, "tau_M"),
        }
        super().__init__(converter, subsystems, connections, pwm=pwm, delay=delay)
