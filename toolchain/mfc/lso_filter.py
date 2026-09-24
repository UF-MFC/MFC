"""
LSO (Least-Squares Optimal) variable-weight 9-point filter for MFC.

Given the particle diameter (d_p), grid spacing (dx/dy/dz), and target
Gaussian filter width (filter_sigma), this module computes the per-pass
filter coefficients via Block Coordinate Descent (BCD) optimization and
returns them ready for insertion into the Fortran namelist.

Filter transfer function (9-point symmetric stencil):
    H(xi) = a0 + 2*a1*cos(xi) + 2*a2*cos(2*xi) + 2*a3*cos(3*xi) + 2*a4*cos(4*xi)

Target (Gaussian in frequency domain):
    G(xi) = exp(-0.5 * sigma_target^2 * xi^2)
    where sigma_target = filter_sigma / dx  (in grid-cell units)

Normalisation constraint (mass conservation):
    H(0) = 1  =>  a0 + 2*(a1 + a2 + a3 + a4) = 1

References:
    Vanden & Belcourt (1995), BC Kim implementation (2026)
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

# ── tuneable constants (must match Fortran side) ──────────────────────────────
LSO_MAX_PASSES: int = 60  # upper bound on filter passes; matches lso_max_passes in Fortran
# A single cascade becomes poorly conditioned for very wide targets; split those
# targets into two equal-variance cascades on the same grid.
PP_SPLIT_CELLS: float = 30.0
PP_STAGE_MAX_CELLS: float = 40.0
_CONV_TOL: float = 1e-3  # default frequency-domain L2 convergence tolerance
_N_XI: int = 600  # quadrature points over [0, pi]
# ─────────────────────────────────────────────────────────────────────────────


def _filter_tf(coeffs: np.ndarray, xi: np.ndarray) -> np.ndarray:
    """Evaluate the 9-point filter transfer function at wavenumbers xi."""
    a0, a1, a2, a3, a4 = coeffs
    return a0 + 2.0 * a1 * np.cos(xi) + 2.0 * a2 * np.cos(2.0 * xi) + 2.0 * a3 * np.cos(3.0 * xi) + 2.0 * a4 * np.cos(4.0 * xi)


def _taylor_init(sigma_per_pass: float) -> np.ndarray:
    """
    Compute Taylor-matched initial coefficients for a single pass.

    Solves the 5x5 linear system that matches H(xi) = exp(-r0^2*xi^2/2) to
    O(xi^8) by equating coefficients at xi = 1, 2, 3, 4 (in normalised units)
    plus the normalisation constraint H(0) = 1.
    """
    r0 = min(sigma_per_pass, 1.38)
    r02 = r0**2

    rows = []
    for m in [1, 2, 3, 4]:
        row = [
            1.0,
            -(m * r02) / 2.0,
            (m * r02) ** 2 / 8.0,
            -((m * r02) ** 3) / 48.0,
            (m * r02) ** 4 / 384.0,
        ]
        rows.append(row)
    rows.append([1.0, 2.0, 2.0, 2.0, 2.0])  # normalisation H(0)=1

    A = np.array(rows, dtype=float)
    b = np.array(
        [np.exp(-0.5 * m**2 * r02) for m in [1, 2, 3, 4]] + [1.0],
        dtype=float,
    )
    try:
        return np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return np.full(5, 0.2)


def design_lso_filter(
    sigma_target: float,
    n_passes: int,
    tol_inner: float = 1e-8,
    n_max: int = 500,
    n_xi: int = _N_XI,
) -> List[Tuple[float, ...]]:
    """
    Design a variable-weight LSO filter with n_passes passes via BCD.

    Minimises sum_xi (H_1 * H_2 * ... * H_K - G)^2  subject to H_k(0) = 1
    for each pass k, iterating until the maximum per-pass weight change
    is below tol_inner.

    Args:
        sigma_target: Target Gaussian width in grid-cell units (filter_sigma/dx).
        n_passes:     Number of filter passes.
        tol_inner:    BCD inner convergence tolerance.
        n_max:        Maximum BCD iterations.
        n_xi:         Quadrature points in [0, pi].

    Returns:
        List of n_passes tuples (a0, a1, a2, a3, a4).
    """
    xi = np.linspace(0.0, np.pi, n_xi)
    G = np.exp(-0.5 * sigma_target**2 * xi**2)
    Phi = np.column_stack(
        [
            np.ones(n_xi),
            2.0 * np.cos(xi),
            2.0 * np.cos(2.0 * xi),
            2.0 * np.cos(3.0 * xi),
            2.0 * np.cos(4.0 * xi),
        ]
    )
    c = np.array([1.0, 2.0, 2.0, 2.0, 2.0])  # H(0)=1 constraint vector

    sigma_per_pass = sigma_target / np.sqrt(max(n_passes, 1))
    betas = [_taylor_init(sigma_per_pass).copy() for _ in range(n_passes)]

    for _ in range(n_max):
        prev = [b.copy() for b in betas]

        for k in range(n_passes):
            H_rest = np.ones(n_xi)
            for j, b in enumerate(betas):
                if j != k:
                    H_rest *= Phi @ b

            Phi_w = H_rest[:, None] * Phi
            A_w = Phi_w.T @ Phi_w + 1e-9 * np.eye(5)
            rhs = Phi_w.T @ G
            A_inv = np.linalg.inv(A_w)
            beta_u = A_inv @ rhs
            Ainv_c = A_inv @ c
            lam = (c @ beta_u - 1.0) / (c @ Ainv_c + 1e-12)
            betas[k] = beta_u - lam * Ainv_c

        delta = max(np.linalg.norm(betas[k] - prev[k]) for k in range(n_passes))
        if delta < tol_inner:
            break

    return [tuple(float(x) for x in b) for b in betas]


def find_min_lso_passes(
    sigma_target: float,
    conv_tol: float = _CONV_TOL,
    max_passes: int = LSO_MAX_PASSES,
    n_xi: int = _N_XI,
) -> Tuple[int, float, List[Tuple[float, ...]]]:
    """
    Find the minimum number of passes to achieve conv_tol.

    The error metric is IC-independent (frequency-domain L2):
        err = sqrt( mean( (H_cum(xi) - G(xi))^2 ) )
    over xi in [0, pi].

    Args:
        sigma_target: Target Gaussian width in grid-cell units.
        conv_tol:     Frequency-domain L2 convergence tolerance.
        max_passes:   Maximum number of passes to try.
        n_xi:         Quadrature points.

    Returns:
        (n_passes, err, coeffs_list)
    """
    xi = np.linspace(0.0, np.pi, n_xi)
    G = np.exp(-0.5 * sigma_target**2 * xi**2)

    err = float("inf")
    coeffs: List[Tuple[float, ...]] = []

    for n_passes in range(1, max_passes + 1):
        coeffs = design_lso_filter(sigma_target, n_passes, n_xi=n_xi)
        H_cum = np.ones(n_xi)
        for ck in coeffs:
            H_cum *= _filter_tf(np.array(ck), xi)
        err = float(np.sqrt(np.mean((H_cum - G) ** 2)))
        if err < conv_tol:
            return n_passes, err, coeffs

    raise ValueError(f"LSO design did not reach L2 tolerance {conv_tol:.3e} in {max_passes} passes; final error is {err:.3e}")


def compute_lso_params(
    d_p: float,
    dx: float,
    dy: float,
    dz: float,
    filter_sigma: float,
    conv_tol: float = _CONV_TOL,
    max_passes: int = LSO_MAX_PASSES,
) -> dict:
    """
    Compute all LSO filter parameters from physical simulation inputs.

    This is the main entry point called by case.py when lso_filter = T.

    Args:
        d_p:          Particle diameter (physical units, same as domain coords).
        dx, dy, dz:   Grid spacing in each direction (physical units).
                      Pass dz = 0.0 for 2D simulations (z-filter is skipped).
        filter_sigma: Target Gaussian filter standard deviation (physical units).
        conv_tol:     Frequency-domain L2 convergence tolerance.
        max_passes:   Upper bound on filter passes; must equal LSO_MAX_PASSES.

    Returns:
        dict with keys:
            lso_n_passes_x/y/z  (int; inactive directions are zero)
            lso_a_x/y/z         (list of n_passes tuples of 5 floats)
    """
    if d_p <= 0.0 or filter_sigma <= 0.0 or dx <= 0.0:
        raise ValueError("d_p, filter_sigma, and dx must be positive")

    directions = [("x", dx)]
    if dy > 0.0:
        directions.append(("y", dy))
    if dz > 0.0:
        directions.append(("z", dz))

    result: dict = {}
    for tag, d in directions:
        sigma_target = filter_sigma / d  # convert physical -> grid-cell units
        n_passes, err, coeffs = find_min_lso_passes(sigma_target, conv_tol=conv_tol, max_passes=max_passes)
        result[f"lso_n_passes_{tag}"] = n_passes
        result[f"lso_a_{tag}"] = coeffs
        n_res = d_p / d
        print(f"  [LSO] {tag}-dir: N_res={n_res:.2f}, sigma_target={sigma_target:.2f} cells, " f"N_passes={n_passes}, L2_err={err:.2e}")

    # Inactive directions use zero passes and dummy weights.
    if dy <= 0.0:
        result["lso_n_passes_y"] = 0
        result["lso_a_y"] = []
    if dz <= 0.0:
        result["lso_n_passes_z"] = 0
        result["lso_a_z"] = []

    return result


def lso_namelist_lines(lso_params: dict, prefix: str = "lso") -> str:
    """
    Format LSO parameters as Fortran namelist lines.

    The weight array lso_a_x has shape (5, LSO_MAX_PASSES) in Fortran
    (column-major). We write only the non-zero entries (up to n_passes per
    direction) to keep the namelist compact; the rest are zero by default.

    Args:
        lso_params: dict returned by compute_lso_params().

    Returns:
        Multi-line string ready to be appended to the &user_inputs namelist.
    """
    lines = []

    for tag in ["x", "y", "z"]:
        n_passes = lso_params.get(f"lso_n_passes_{tag}", 0)
        coeffs = lso_params.get(f"lso_a_{tag}", [])

        lines.append(f"{prefix}_n_passes_{tag} = {n_passes}")

        if n_passes > 0 and coeffs:
            # Build flat array in Fortran column-major order:
            # lso_a_x(coeff_idx, pass_idx) with coeff_idx varying fastest
            flat = []
            for i_pass in range(LSO_MAX_PASSES):
                for j_coeff in range(5):
                    if i_pass < len(coeffs):
                        flat.append(f"{coeffs[i_pass][j_coeff]:.17e}")
                    else:
                        flat.append("0.0d0")
            lines.append(f"{prefix}_a_{tag} = {' '.join(flat)}")

    return "\n".join(lines) + "\n"
