"""
Utilities for multi-three-phase machine drives.

The machine has n three-phase winding sets. The phase a of the set k is located at the
angle `theta_sets[k]` with respect to the stationary α axis. The space vector of the set
k in its own coordinates is `abc2complex(x_abc_k)`, and in the common stationary
coordinates `exp(1j*theta_sets[k])*abc2complex(x_abc_k)`.

The decoupled multi-stator (DMS) transformation maps the space vectors of the sets in
the common stationary coordinates to the common mode (CM) and n - 1 differential modes
(DMs), see :func:`dms_matrix`. For two sets, `x_cm = (x_1 + x_2)/2` and `x_dm = (x_1 -
x_2)/2`. Including the set angles, the modes are `D @ x_sets`, where `D =
T*exp(1j*theta_sets)` and `x_sets` are the set vectors in their own coordinates, and the
inverse is `n*D.conj().T`.

The plotting function and the base values of motulator are re-exported for
convenience.

"""

from collections.abc import Callable
from dataclasses import replace
from functools import partial
from math import pi
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from motulator.common.model import SimulationResults
from motulator.common.utils import set_style
from motulator.drive.model import (
    InductionMachineInvGammaPars,
    InductionMachinePars,
    SaturatedSynchronousMachinePars,
    SpatialSaturatedSynchronousMachinePars,
    SynchronousMachinePars,
)
from motulator.drive.utils import BaseValues, NominalValues, plot

__all__ = [
    "BaseValues",
    "MachinePars",
    "NominalValues",
    "default_set_angles",
    "dms_matrix",
    "plot",
    "plot_sets",
    "to_parallel",
]

# %%
type MachinePars = (
    SynchronousMachinePars
    | SaturatedSynchronousMachinePars
    | SpatialSaturatedSynchronousMachinePars
    | InductionMachinePars
    | InductionMachineInvGammaPars
)


def default_set_angles(n_sets: int) -> np.ndarray:
    """
    Angles of the winding sets of an asymmetrical multi-three-phase machine.

    The sets are displaced by `pi/(3*n_sets)`, e.g., 30° for a dual-three-phase machine
    and 15° for a quadruple-three-phase machine.

    Parameters
    ----------
    n_sets : int
        Number of three-phase sets.

    Returns
    -------
    ndarray, shape (n_sets,)
        Angles (rad) of the phases a of the sets.

    Examples
    --------
    >>> from motulator_multiphase.utils import default_set_angles
    >>> np.degrees(default_set_angles(2))
    array([ 0., 30.])

    """
    return np.arange(n_sets) * pi / (3 * n_sets)


def dms_matrix(n_sets: int) -> np.ndarray:
    """
    Decoupled multi-stator (DMS) transformation matrix.

    The transformation maps the space vectors of the sets in common coordinates to the
    common mode and the differential modes [#Rub2021]_. The first row gives the common
    mode as the mean of the sets. The other rows give the differential modes, `x_dm_u ∝
    (n - u)*x_u - sum(x_{u+1}, ..., x_n)`. The rows are orthogonal and scaled such that
    the inverse is `n_sets*T.T`.

    Parameters
    ----------
    n_sets : int
        Number of three-phase sets.

    Returns
    -------
    ndarray, shape (n_sets, n_sets)
        Transformation matrix.

    References
    ----------
    .. [#Rub2021] Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation
       for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power
       Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083

    Examples
    --------
    >>> from motulator_multiphase.utils import dms_matrix
    >>> dms_matrix(2)
    array([[ 0.5,  0.5],
           [ 0.5, -0.5]])
    >>> T = dms_matrix(4)
    >>> bool(np.allclose(4 * T.T @ T, np.eye(4)))
    True

    """
    n = n_sets
    T = np.zeros((n, n))
    T[0] = 1 / n
    for u in range(1, n):
        row = np.zeros(n)
        row[u - 1] = n - u
        row[u:] = -1
        T[u] = row / np.sqrt(n * row @ row)
    return T


# %%
# The scaled maps are module-level functions with bound arguments, so that the models
# and control systems can be pickled (e.g., for parallel simulations)
def _scaled_jacobian(
    jacobian: Callable[[Any], Any], k_in: float, k_out: float, x: Any
) -> np.ndarray:
    return k_out * k_in * np.asarray(jacobian(k_in * x))


def _scaled_magnetic_map(
    fcn: Callable[[Any, Any], tuple[Any, Any]],
    n: int,
    psi_s_dq: Any,
    exp_j_theta_m: Any,
) -> tuple[Any, Any]:
    i_s_dq, tau_m = fcn(psi_s_dq, exp_j_theta_m)
    return n * i_s_dq, n * tau_m


def _scaled_inductance(L_s: Callable[[float], float], n: int, psi: float) -> float:
    return L_s(psi) / n


class _ScaledMap:
    """Map `x -> k_out*fcn(k_in*x)`, with the Jacobian of `fcn` if it has one."""

    def __init__(self, fcn: Callable[[Any], Any], k_in: float, k_out: float) -> None:
        self.fcn = fcn
        self.k_in = k_in
        self.k_out = k_out
        jacobian = getattr(fcn, "jacobian", None)
        self.jacobian: Callable[[Any], np.ndarray] | None = (
            None
            if jacobian is None
            else partial(_scaled_jacobian, jacobian, k_in, k_out)
        )

    def __call__(self, x: Any) -> Any:
        return self.k_out * self.fcn(self.k_in * x)


def to_parallel[P: MachinePars](par: P, n_sets: int) -> P:
    """
    Parameters of the equivalent machine with the sets connected in parallel.

    The given parameters describe one three-phase set, when all sets carry equal
    currents (i.e., the common mode only). When the sets are connected in parallel,
    the voltages and the flux linkages remain unchanged, while the currents are
    multiplied by `n_sets`. Hence, the resistances and inductances are divided by
    `n_sets`, and the core-loss conductance is multiplied by it. The equivalent machine
    produces the total torque of the machine. The maps of saturated machines are scaled
    accordingly, keeping their analytic Jacobians (e.g., of GradNet flux maps).

    Parameters
    ----------
    par : MachinePars
        Parameters of one set (balanced operation of all sets).
    n_sets : int
        Number of three-phase sets.

    Returns
    -------
    MachinePars
        Parameters of the equivalent three-phase machine, of the same type as `par`.

    Examples
    --------
    >>> from motulator_multiphase import model, utils
    >>> par = model.SynchronousMachinePars(
    ...     n_p=3, R_s=3.6, L_d=0.036, L_q=0.051, psi_f=0.545
    ... )
    >>> utils.to_parallel(par, 2)
    SynchronousMachinePars(n_p=3, R_s=1.8, L_d=0.018, L_q=0.0255, psi_f=0.545, G_c=0.0)

    """
    n = n_sets
    R_s, G_c = par.R_s / n, n * par.G_c
    if isinstance(par, SynchronousMachinePars):
        return replace(par, R_s=R_s, L_d=par.L_d / n, L_q=par.L_q / n, G_c=G_c)
    if isinstance(par, SaturatedSynchronousMachinePars):
        i_fcn, psi_fcn = par.i_s_dq_fcn, par.psi_s_dq_fcn
        return replace(
            par,
            R_s=R_s,
            i_s_dq_fcn=None if i_fcn is None else _ScaledMap(i_fcn, 1, n),
            psi_s_dq_fcn=None if psi_fcn is None else _ScaledMap(psi_fcn, 1 / n, 1),
            G_c=G_c,
        )
    if isinstance(par, SpatialSaturatedSynchronousMachinePars):
        fcn = partial(_scaled_magnetic_map, par.magnetic_map_fcn, n)
        return replace(par, R_s=R_s, magnetic_map_fcn=fcn, G_c=G_c)
    if isinstance(par, InductionMachinePars):
        L_s = par.L_s
        L_s_n: float | Callable[[float], float] = (
            partial(_scaled_inductance, L_s, n) if callable(L_s) else L_s / n
        )
        return replace(
            par, R_s=R_s, R_r=par.R_r / n, L_ell=par.L_ell / n, L_s=L_s_n, G_c=G_c
        )
    return replace(
        par, R_s=R_s, R_R=par.R_R / n, L_sgm=par.L_sgm / n, L_M=par.L_M / n, G_c=G_c
    )


# %%
def plot_sets(
    res: SimulationResults,
    t_lims: tuple[float, float] | None = None,
    latex: bool = False,
    save_path: str | Path | None = None,
    **savefig_kwargs,
) -> None:
    """
    Plot the set and differential-mode currents.

    The first subplot shows the current magnitudes of the sets. The second subplot shows
    the DM currents in the controller coordinates (d solid, q dash-dotted) and their
    references (black dashed).

    Parameters
    ----------
    res : SimulationResults
        Simulation results of a multi-three-phase drive.
    t_lims : tuple[float, float], optional
        Time axis limits. If None, uses full time range.
    latex : bool, optional
        Use LaTeX fonts for the labels, defaults to False.
    save_path : str | Path, optional
        Path to save the figure. If None, the figure is not saved.
    **savefig_kwargs
        Additional keyword arguments passed to `plt.savefig()`.

    """
    set_style(latex)
    mdl, ctrl = res.mdl, res.ctrl
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    for k, i_k in enumerate(np.abs(mdl.machine.i_c_ab)):
        ax1.plot(mdl.t, i_k, label=rf"$|\mathbf{{i}}_{{{k + 1}}}|$")
    ax1.legend(loc="best", ncols=2)
    for u, (i_dm, i_dm_ref) in enumerate(
        zip(ctrl.dm.i.T, ctrl.dm.i_ref.T, strict=True)
    ):
        c, sub = f"C{u}", rf"\mathrm{{dm}}{u + 1}"
        ax2.plot(ctrl.t, i_dm.real, c=c, ds="steps-post", label=rf"$i_{{{sub},d}}$")
        ax2.plot(
            ctrl.t, i_dm.imag, "-.", c=c, ds="steps-post", label=rf"$i_{{{sub},q}}$"
        )
        ax2.plot(ctrl.t, i_dm_ref.real, "k--", lw=0.5, ds="steps-post")
        ax2.plot(ctrl.t, i_dm_ref.imag, "k--", lw=0.5, ds="steps-post")
    if mdl.machine.psi_dm_ab.shape[0] > 0:  # No DMs with one set
        ax2.legend(loc="best", ncols=2)
    t_lims = (0, mdl.t[-1]) if t_lims is None else t_lims
    ax1.set_xlim(t_lims)
    ax1.set_ylabel("Set current (A)")
    ax2.set_ylabel("DM current (A)")
    fig.align_ylabels()
    ax2.set_xlabel("Time (s)")
    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, **savefig_kwargs)
    plt.show()
