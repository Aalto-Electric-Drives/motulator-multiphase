"""Helper functions copied from the private modules of motulator."""

from collections.abc import Callable
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler


def wrap(theta: float) -> Any:
    """Limit the angle into the range [-pi, pi)."""
    return np.mod(theta + np.pi, 2 * np.pi) - np.pi


def get_value(u: Any | Callable[[Any], Any] | None, x: Any) -> Any:
    """Return `u(x)` if `u` is callable, otherwise `u`."""
    if u is None:
        raise ValueError("Input value cannot be None")
    return u(x) if callable(u) else u


def empty_array() -> np.ndarray:
    """Return an empty array."""
    return np.array([])


def set_screen_style() -> None:
    """Configure matplotlib for screen viewing."""
    plt.style.use("default")
    # plt.style.use("dark_background")
    plt.rcParams.update(
        {
            "axes.prop_cycle": cycler("color", ["b", "r", "m", "c"]),
            "lines.linewidth": 1,
            "axes.grid": True,
            "text.usetex": False,
            "axes3d.grid": True,
            "figure.constrained_layout.use": True,
        }
    )


def set_latex_style() -> None:
    """Configure matplotlib for LaTeX documents."""
    plt.style.use("default")
    plt.rcParams.update(
        {
            "axes.prop_cycle": cycler("color", ["b", "r", "m", "c"]),
            "lines.linewidth": 1,
            "axes.grid": True,
            # "axes.autolimit_mode": "round_numbers",
            "text.usetex": True,
            "text.latex.preamble": (r"\usepackage{newtxtext}\usepackage{newtxmath}"),
            "font.family": "serif",
            "font.serif": ["Liberation Serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 8,
            "axes.labelsize": 8,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes3d.grid": True,
            "figure.figsize": [3.5, 2.5],
            # "figure.constrained_layout.use": True,
            # "savefig.bbox": "tight",
            # "savefig.pad_inches": 0.01,
            # Disable automatic adjustment
            "figure.constrained_layout.use": False,
            "figure.subplot.left": 0.18,  # Fixed margins
            "figure.subplot.right": 0.95,
            "figure.subplot.bottom": 0.18,
            "figure.subplot.top": 0.95,
        }
    )


def set_style(latex: bool) -> None:
    """Configure matplotlib for LaTeX documents or for screen viewing."""
    if latex:
        set_latex_style()
    else:
        set_screen_style()
