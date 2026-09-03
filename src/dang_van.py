"""Educational implementation of a macroscopic Dang Van fatigue check.

Stress components use the order [sxx, syy, szz, sxy, syz, sxz] in MPa.
The implementation deliberately exposes every modelling step so that it can be
reviewed and adapted to a validated industrial workflow.
"""

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass(frozen=True)
class DangVanMaterial:
    """Material limits used by the Dang Van line tau + a*p <= b."""

    reversed_bending_limit: float
    reversed_torsion_limit: float

    def __post_init__(self) -> None:
        if self.reversed_bending_limit <= 0 or self.reversed_torsion_limit <= 0:
            raise ValueError("Fatigue limits must be strictly positive.")

    @property
    def a(self) -> float:
        return 3.0 * self.reversed_torsion_limit / self.reversed_bending_limit - 1.5

    @property
    def b(self) -> float:
        return self.reversed_torsion_limit


@dataclass(frozen=True)
class DangVanResult:
    equivalent_stress: float
    safety_factor: float
    critical_normal: np.ndarray
    critical_direction: np.ndarray
    critical_time_index: int


def voigt_to_tensor(stress_history: np.ndarray) -> np.ndarray:
    """Convert an (N, 6) Voigt history into N symmetric 3x3 tensors."""
    stress_history = np.asarray(stress_history, dtype=float)
    if stress_history.ndim != 2 or stress_history.shape[1] != 6:
        raise ValueError("stress_history must have shape (N, 6).")
    tensors = np.zeros((len(stress_history), 3, 3), dtype=float)
    tensors[:, 0, 0], tensors[:, 1, 1], tensors[:, 2, 2] = stress_history[:, :3].T
    tensors[:, 0, 1] = tensors[:, 1, 0] = stress_history[:, 3]
    tensors[:, 1, 2] = tensors[:, 2, 1] = stress_history[:, 4]
    tensors[:, 0, 2] = tensors[:, 2, 0] = stress_history[:, 5]
    return tensors


def hydrostatic_stress(tensors: np.ndarray) -> np.ndarray:
    """Return p(t) = trace(sigma(t))/3 in MPa."""
    return np.trace(tensors, axis1=1, axis2=2) / 3.0


def _orientations(n_theta: int, n_phi: int, n_psi: int):
    """Generate plane normals n and in-plane shear directions m."""
    for theta, phi in product(
        np.linspace(0.0, np.pi / 2.0, n_theta, endpoint=True),
        np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False),
    ):
        n = np.array([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ])
        e1 = np.array([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)])
        e2 = np.array([-np.sin(phi), np.cos(phi), 0.0])
        for psi in np.linspace(0.0, np.pi, n_psi, endpoint=False):
            yield n, np.cos(psi) * e1 + np.sin(psi) * e2


def evaluate_dang_van(
    stress_history: np.ndarray,
    material: DangVanMaterial,
    n_theta: int = 10,
    n_phi: int = 20,
    n_psi: int = 18,
) -> DangVanResult:
    """Search the most damaging plane/direction over a periodic stress history.

    The alternating resolved shear is obtained by removing the midrange of the
    resolved shear history. The maximum of |tau_alt(t)| + a*p(t) is then sought
    over time and sampled orientations.
    """
    if min(n_theta, n_phi, n_psi) < 2:
        raise ValueError("Orientation discretisation values must be >= 2.")

    tensors = voigt_to_tensor(stress_history)
    pressure = hydrostatic_stress(tensors)
    best_value = -np.inf
    best_n = best_m = None
    best_time = 0

    for normal, direction in _orientations(n_theta, n_phi, n_psi):
        resolved_shear = np.einsum("i,tij,j->t", direction, tensors, normal)
        midrange = 0.5 * (resolved_shear.max() + resolved_shear.min())
        alternating_shear = np.abs(resolved_shear - midrange)
        indicator = alternating_shear + material.a * pressure
        time_index = int(np.argmax(indicator))
        if indicator[time_index] > best_value:
            best_value = float(indicator[time_index])
            best_n = normal.copy()
            best_m = direction.copy()
            best_time = time_index

    safety_factor = np.inf if best_value <= 0 else material.b / best_value
    return DangVanResult(best_value, float(safety_factor), best_n, best_m, best_time)


def synthetic_non_proportional_history(n_points: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """Create a generic out-of-phase tension/torsion cycle in MPa."""
    phase = np.linspace(0.0, 2.0 * np.pi, n_points)
    stress = np.zeros((n_points, 6))
    stress[:, 0] = 80.0 + 140.0 * np.sin(phase)
    stress[:, 3] = 75.0 * np.sin(phase + np.pi / 2.0)
    return phase, stress
