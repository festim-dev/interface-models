"""Closed-form steady state of Model 2 on a metal/liquid two-slab problem.

The metal occupies [0, L_m] with an atomic species H, the liquid occupies
[L_m, L_m + L_s] with the molecular carrier H2, and the interface at x = L_m
carries the single recombination channel

    2 H(m) <-> H2(s),   w = k_plus * c_m|G**2 - k_minus * c_s|G

with interface conditions -D_m dn c_m = 2w and D_s dn c_s = w (Eq. model2_rate
and the two conditions following it in the paper).

At steady state both profiles are linear and the same atomic flux phi = 2w
crosses the metal while the molecular flux w crosses the liquid, so

    c_m|G = c_0 - 2 w R_m,   R_m = L_m / D_m
    c_s|G = c_L +   w R_s,   R_s = L_s / D_s

Substituting into the rate leaves a quadratic in w, which is solved exactly.
Unlike Model 1 the result is not a sum of resistances in series: the
interfacial resistance depends on the loading it is evaluated at.

`solution` is the kinetic problem at finite k_plus; `lte_solution` is the
k_plus -> infinity limit at fixed K = k_plus/k_minus, i.e. the Sieverts/Henry
closure c_s|G = K c_m|G**2 (Eq. lte_sh). The sweep of the paper checks the
first against FESTIM and its convergence to the second.
"""

import numpy as np


def _physical_root(a, b, c):
    """Root of a w**2 + b w + c = 0 that stays finite as a -> 0.

    b is negative for every case here (see below), so -b + sqrt(disc) never
    cancels and the citardauq form is both stable and continuous in the limit
    a -> 0, where it returns -c/b.
    """
    disc = b**2 - 4 * a * c
    if disc < 0:
        raise ValueError(f"no real steady state: discriminant {disc:.3e} < 0")
    return 2 * c / (-b + np.sqrt(disc))


def solution(c_0, c_L, D_m, D_s, L_m, L_s, k_plus, k_minus):
    """Steady state of the kinetic channel.

    Args:
        c_0: atomic concentration imposed on the outer metal face
        c_L: molecular concentration imposed on the outer liquid face
        D_m, D_s: diffusivities in the metal and in the liquid
        L_m, L_s: slab thicknesses
        k_plus, k_minus: forward and reverse constants of the channel

    Returns:
        (c_m_int, c_s_int, w) with w the molecular rate; the atomic flux
        leaving the metal is 2w.
    """
    R_m = L_m / D_m
    R_s = L_s / D_s

    a = 4 * k_plus * R_m**2
    b = -(1 + 4 * k_plus * c_0 * R_m + k_minus * R_s)
    c = k_plus * c_0**2 - k_minus * c_L

    w = _physical_root(a, b, c)

    return c_0 - 2 * w * R_m, c_L + w * R_s, w


def lte_solution(c_0, c_L, D_m, D_s, L_m, L_s, K):
    """Steady state of the Sieverts/Henry closure, c_s|G = K c_m|G**2.

    This is the k_plus -> infinity limit of `solution` at fixed
    K = k_plus/k_minus, which detailed balance fixes at K_H/K_S**2
    (Eq. detailed_balance_2).
    """
    R_m = L_m / D_m
    R_s = L_s / D_s

    a = 4 * K * R_m**2
    b = -(4 * K * c_0 * R_m + R_s)
    c = K * c_0**2 - c_L

    w = _physical_root(a, b, c)

    return c_0 - 2 * w * R_m, c_L + w * R_s, w
