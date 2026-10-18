"""Helper functions copied from the private modules of motulator."""

from collections.abc import Callable
from typing import Any

import numpy as np


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
