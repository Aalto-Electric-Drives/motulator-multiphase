"""
Tests for the multi-three-phase machine drives.

The machine models are checked against the multi-stator (MS) models, in which the flux
linkage of each set is computed directly from the currents of all sets and the rotor,
and the drive with one set against the three-phase drive. The closed-loop tests check
speed tracking, torque balance, and current sharing with loose tolerances.

"""

import pickle
from math import pi
from typing import Any

import motulator.drive.control.im as im_control
import motulator.drive.control.sm as sm_control
import numpy as np
import pytest
from motulator.drive import model as drive_model
from scipy.integrate import solve_ivp

from motulator_multiphase import control, model, utils

T_STOP = 0.6
T_LOAD = 0.3
T_SHARE = 0.45
SM_PAR = model.SynchronousMachinePars(n_p=3, R_s=3.6, L_d=0.036, L_q=0.051, psi_f=0.545)
IM_PAR = model.InductionMachineInvGammaPars(
    n_p=2, R_s=3.7, R_R=2.1, L_sgm=0.021, L_M=0.224
)


@pytest.mark.parametrize("n_sets", [1, 2, 3, 4, 5])
def test_dms_matrix(n_sets: int) -> None:
    """The CM is the mean of the sets and the inverse is n*T.T."""
    T = utils.dms_matrix(n_sets)
    assert T[0] == pytest.approx(np.full(n_sets, 1 / n_sets))
    assert n_sets * T.T @ T == pytest.approx(np.eye(n_sets))


class _SaturatingFluxMap:
    """Analytic saturating flux map psi = L0*i/(1 + a*|i|**2) with its Jacobian."""

    L0, a = 0.01, 0.002

    def __call__(self, i_s_dq: Any) -> Any:
        return self.L0 * i_s_dq / (1 + self.a * np.abs(i_s_dq) ** 2)

    def jacobian(self, i_s_dq: Any) -> np.ndarray:
        i = np.asarray(i_s_dq)
        x, y, g = i.real, i.imag, 1 + self.a * np.abs(i) ** 2
        k = 2 * self.a * self.L0 / g**2
        J = np.array(
            [
                [self.L0 / g - k * x * x, -k * x * y],
                [-k * x * y, self.L0 / g - k * y * y],
            ]
        )
        return np.moveaxis(J, (0, 1), (-2, -1))


def test_parallel_flux_map_keeps_the_jacobian() -> None:
    """The analytic Jacobian of a flux map is scaled, not replaced by differences."""
    par = model.SaturatedSynchronousMachinePars(
        n_p=2, R_s=0.63, psi_s_dq_fcn=_SaturatingFluxMap()
    )
    i_s_dq = np.array([1 + 2j, -3 + 0.5j])
    L_s = utils.to_parallel(par, 4).incr_ind_mat(i_s_dq)
    assert L_s == pytest.approx(par.incr_ind_mat(i_s_dq / 4) / 4, rel=1e-12)


@pytest.mark.parametrize("n_sets", [2, 3])
def test_machine_matches_multi_stator_model(n_sets: int) -> None:
    """Set currents and torque agree with the MS model of a linear PMSM."""
    par, L_ls = SM_PAR, 0.005
    machine = model.MultiThreePhaseMachine(par, L_ls=L_ls, n_sets=n_sets)
    # MS model: psi_k = L_ls*i_k + M*sum(i_z) + psi_f, rotor aligned with the α axis
    M_d, M_q = (par.L_d - L_ls) / n_sets, (par.L_q - L_ls) / n_sets
    rng = np.random.default_rng(1)
    i_s = rng.normal(size=n_sets) + 1j * rng.normal(size=n_sets)
    i_sum = np.sum(i_s)
    psi_s = L_ls * i_s + M_d * i_sum.real + 1j * M_q * i_sum.imag + par.psi_f
    # Set the states: CM flux in rotor coordinates and DM fluxes in stator coordinates
    T = utils.dms_matrix(n_sets)
    psi_modes = T @ psi_s
    machine.set_states([psi_modes[0], 0j, *psi_modes[1:]], 0)
    machine.set_outputs(0.0)

    i_s_model = np.exp(1j * machine.theta_sets) * machine.out.i_c_ab
    tau_M = 1.5 * par.n_p * np.sum(np.imag(i_s * np.conj(psi_s)))
    assert i_s_model == pytest.approx(i_s)
    assert machine.out.tau_M == pytest.approx(tau_M)


@pytest.mark.parametrize("n_sets", [2, 3])
def test_induction_machine_matches_multi_stator_model(n_sets: int) -> None:
    """Set currents and torque agree with the MS T-model of an induction machine."""
    R_s, L_ls, L_m, L_lr, R_r, n_p, w_M = 3.7, 0.008, 0.1, 0.006, 1.0, 2, 100.0
    rng = np.random.default_rng(0)
    # Unequal set voltages, which excite also the DMs (common stationary coordinates)
    U = (
        200
        * (1 + 0.3 * rng.normal(size=n_sets))
        * np.exp(0.4j * rng.normal(size=n_sets))
    )

    def u_s(t: float) -> np.ndarray:
        return U * np.exp(2j * pi * 40 * t)

    t_eval = np.linspace(0, 0.1, 101)
    tol = {"t_eval": t_eval, "rtol": 1e-9, "atol": 1e-12}

    # MS T-model: psi_k = L_ls*i_k + L_m*(sum(i_z) + i_r), psi_r = L_lr*i_r + L_m*(...)
    L = L_m * np.ones((n_sets + 1, n_sets + 1)) + np.diag(n_sets * [L_ls] + [L_lr])

    def rhs_ms(t: float, psi: np.ndarray) -> np.ndarray:
        i = np.linalg.solve(L, psi)
        d_psi_s = u_s(t) - R_s * i[:n_sets]
        d_psi_r = -R_r * i[n_sets] + 1j * n_p * w_M * psi[n_sets]
        return np.append(d_psi_s, d_psi_r)

    psi0 = np.zeros(n_sets + 1, dtype=complex)
    sol_ms = solve_ivp(rhs_ms, (0, t_eval[-1]), psi0, **tol)
    i_ms = np.linalg.solve(L, sol_ms.y)[:n_sets]
    tau_ms = 1.5 * n_p * np.sum(np.imag(i_ms * np.conj(sol_ms.y[:n_sets])), axis=0)

    # Balanced set: T-model with L_ls, n*L_m, n*L_lr, and n*R_r, mapped to the Γ model
    L_M, L_s = n_sets * L_m, L_ls + n_sets * L_m
    k = L_s / L_M
    par = model.InductionMachinePars(
        n_p=n_p,
        R_s=R_s,
        R_r=k**2 * n_sets * R_r,
        L_ell=k**2 * n_sets * L_lr + k * L_ls,
        L_s=L_s,
    )
    machine = model.MultiThreePhaseMachine(par, L_ls=L_ls, n_sets=n_sets)
    rot = np.exp(1j * machine.theta_sets)
    machine.inp.w_M = w_M

    def outputs(t: float, y: np.ndarray) -> None:
        machine.set_states(list(y), 0)
        machine.inp.u_c_ab = u_s(t) / rot
        machine.set_outputs(t)

    def rhs(t: float, y: np.ndarray) -> list[complex]:
        outputs(t, y)
        return machine.rhs(t)

    y0 = np.zeros(len(vars(machine.state)), dtype=complex)
    sol = solve_ivp(rhs, (0, t_eval[-1]), y0, **tol)
    i_s, tau_M = [], []
    for t, y in zip(sol.t, sol.y.T, strict=True):
        outputs(t, y)
        i_s.append(rot * machine.out.i_c_ab)
        tau_M.append(machine.out.tau_M)

    assert np.array(i_s).T == pytest.approx(i_ms, abs=1e-5 * np.max(np.abs(i_ms)))
    assert np.array(tau_M) == pytest.approx(tau_ms, abs=1e-5 * np.max(np.abs(tau_ms)))


def test_one_set_equals_three_phase_drive() -> None:
    """With one set, the drive and its control equal the three-phase ones."""

    def simulate(mdl: Any, ctrl: Any) -> Any:
        mdl.mechanics.set_external_rotor_speed(lambda t: min(t / 0.1, 1) * 100.0)
        ctrl.set_torque_ref(lambda t: (t > 0.05) * 10.0)
        return model.Simulation(mdl, ctrl, show_progress=False).simulate(0.2)

    def vector_ctrl() -> sm_control.CurrentVectorController:
        cfg = sm_control.CurrentVectorControllerCfg(i_s_max=6.5, sensorless=False)
        return sm_control.CurrentVectorController(SM_PAR, cfg)

    res = simulate(
        drive_model.Drive(
            drive_model.SynchronousMachine(SM_PAR),
            model.ExternalRotorSpeed(),
            drive_model.VoltageSourceConverter(u_dc=540),
        ),
        sm_control.VectorControlSystem(vector_ctrl()),
    )
    res_1 = simulate(
        model.MultiThreePhaseDrive(
            model.MultiThreePhaseMachine(SM_PAR, L_ls=0.005, n_sets=1),
            model.ExternalRotorSpeed(),
            model.MultiThreePhaseConverter(u_dc=540, n_sets=1),
        ),
        control.MultiThreePhaseControlSystem(
            vector_ctrl(), control.DifferentialModeController(0.005, n_sets=1)
        ),
    )

    assert res_1.mdl.t == pytest.approx(res.mdl.t)
    assert res_1.mdl.machine.i_c_ab[0] == pytest.approx(res.mdl.machine.i_s_ab)
    assert res_1.mdl.machine.tau_M == pytest.approx(res.mdl.machine.tau_M, abs=1e-6)


def make_case(
    ctrl_type: str, n_sets: int, pwm: bool, discrete: bool = False
) -> tuple[model.MultiThreePhaseDrive, control.MultiThreePhaseControlSystem]:
    """Create the drive model and control system."""
    if ctrl_type == "sm_fvc":
        cfg = sm_control.FluxVectorControllerCfg(i_s_max=n_sets * 6.5)
        est_par = utils.to_parallel(SM_PAR, n_sets)
        vector_ctrl = sm_control.FluxVectorController(est_par, cfg)
    elif ctrl_type == "sm_cvc":
        cfg = sm_control.CurrentVectorControllerCfg(
            i_s_max=n_sets * 6.5, sensorless=False, discrete=discrete
        )
        est_par = utils.to_parallel(SM_PAR, n_sets)
        vector_ctrl = sm_control.CurrentVectorController(est_par, cfg)
    else:
        cfg = im_control.FluxVectorControllerCfg(psi_s_nom=1.04, i_s_max=n_sets * 10.6)
        est_par = utils.to_parallel(IM_PAR, n_sets)
        vector_ctrl = im_control.FluxVectorController(est_par, cfg)
    par, L_ls = (IM_PAR, 0.008) if ctrl_type == "im_fvc" else (SM_PAR, 0.005)
    J = n_sets * 0.015
    mdl = model.MultiThreePhaseDrive(
        model.MultiThreePhaseMachine(par, L_ls=L_ls, n_sets=n_sets),
        model.MechanicalSystem(J=J),
        model.MultiThreePhaseConverter(u_dc=540, n_sets=n_sets),
        pwm=pwm,
    )
    ctrl = control.MultiThreePhaseControlSystem(
        vector_ctrl,
        control.DifferentialModeController(
            L_ls, n_sets=n_sets, T_s=cfg.T_s if discrete else None
        ),
        sm_control.SpeedController(J=J, alpha_s=2 * pi * 4),
    )
    return mdl, ctrl


CASES = {
    "sm_fvc_sensorless": ("sm_fvc", 2, False),
    "sm_cvc_sensored": ("sm_cvc", 2, False),
    "sm_cvc_pwm": ("sm_cvc", 2, True),
    "sm_cvc_discrete": ("sm_cvc", 3, False, True),
    "sm_fvc_four_sets": ("sm_fvc", 4, False),
    "im_fvc_sensorless": ("im_fvc", 2, False),
}


@pytest.mark.parametrize("case", CASES)
def test_speed_control_and_sharing(case: str) -> None:
    """Speed follows the reference, torque balances the load, currents are shared."""
    mdl, ctrl = make_case(*CASES[case])
    n_sets = mdl.machine.n_sets
    w_M_ref, tau_L = 50.0, n_sets * 7.0
    k_share = np.ones(n_sets)
    k_share[:2] += [0.5, -0.5]

    ctrl.set_speed_ref(lambda t: (t > 0.05) * w_M_ref)
    mdl.mechanics.set_external_load_torque(lambda t: (t > T_LOAD) * tau_L)
    ctrl.set_current_sharing(lambda t: k_share if t > T_SHARE else np.ones(n_sets))
    res = model.Simulation(mdl, ctrl, show_progress=False).simulate(t_stop=T_STOP)

    # Averages over the last 50 ms and before the sharing step
    end = res.mdl.t > T_STOP - 0.05
    before = (res.mdl.t > T_SHARE - 0.05) & (res.mdl.t < T_SHARE)
    ctrl_end = res.ctrl.t > T_STOP - 0.05
    w_M = np.mean(res.mdl.machine.w_M[end])
    tau_M = np.mean(res.mdl.machine.tau_M[end])
    i_sets_before = np.mean(np.abs(res.mdl.machine.i_c_ab[:, before]), axis=1)
    i_sets = np.mean(np.abs(res.mdl.machine.i_c_ab[:, end]), axis=1)
    i_dm = np.mean(res.ctrl.dm.i[ctrl_end], axis=0)
    i_dm_ref = np.mean(res.ctrl.dm.i_ref[ctrl_end], axis=0)

    assert w_M == pytest.approx(w_M_ref, rel=0.02)
    assert tau_M == pytest.approx(tau_L, abs=0.05 * n_sets * 14)
    assert np.mean(res.ctrl.fbk.w_M[ctrl_end]) == pytest.approx(w_M, rel=0.02)
    # Equal sharing first, then the set currents follow the sharing coefficients
    assert i_sets_before == pytest.approx(np.mean(i_sets_before), rel=0.02)
    assert i_sets / np.mean(i_sets) == pytest.approx(k_share, abs=0.03)
    assert i_dm == pytest.approx(i_dm_ref, abs=0.02 * np.max(i_sets))


def test_models_and_control_systems_can_be_pickled() -> None:
    """Models and control systems can be pickled, e.g., for parallel simulations."""
    for args in CASES.values():
        pickle.dumps(make_case(*args))


@pytest.mark.filterwarnings("ignore:overflow:RuntimeWarning")
def test_numerical_failure_is_flagged() -> None:
    """A failure in the first sampling period is flagged and gives empty results."""
    mdl, ctrl = make_case("sm_fvc", 2, False)
    ctrl.set_speed_ref(0.0)
    mdl.mechanics.set_external_load_torque(1e300)
    sim = model.Simulation(mdl, ctrl, show_progress=False)
    with pytest.warns(RuntimeWarning, match="Simulation stopped"):
        res = sim.simulate(t_stop=0.05)

    assert not res.success
    assert res.mdl.machine.i_c_ab.shape == (2, 0)
