"""Closed-form steady state of Model 3 on a metal/liquid two-slab problem.

The metal occupies [0, L_m] with the atomic species H, the liquid occupies
[L_m, L_m + L_s] with two carriers, H2 and HF, and the interface at x = L_m
carries the two competing channels of Sec. "Model 3: competing channels":

    R:  2 H(m) <-> H2(s),   w_rec = kr_plus * c**2   - kr_minus * c_H2|G
    F:    H(m) <-> HF(s),   w_F   = kf_plus * a_F * c - kf_minus * c_HF|G

with c = c_m|G, a_F the effective fluoride activity folded into the forward
constant as in Eq. (model3_rate), and the atomic flux leaving the metal shared
between the two, -D_m dn c_m = 2 w_rec + w_F (Eq. model3_total_flux). Each
carrier is fed by the single channel that produces it.

At steady state all three profiles are linear, so with R_m = L_m/D_m,
R_2 = L_s/D_H2 and R_F = L_s/D_HF,

    c        = c_0  - (2 w_rec + w_F) R_m
    c_H2|G   = c_L2 +      w_rec      R_2
    c_HF|G   = c_LF +      w_F        R_F

Eliminating the two salt-side traces makes each channel rate an explicit
function of c alone,

    w_rec = A c**2 - a0,   A  = kr_plus / (1 + kr_minus R_2)
    w_F   = B1 c   - b0,   B1 = kf_plus a_F / (1 + kf_minus R_F)

with a0 = kr_minus c_L2 / (1 + kr_minus R_2) and b0 likewise, and the metal
balance then leaves a quadratic in c which is solved exactly.

Two consequences are used by the verification script.

First, the downstream resistance does not change the *form* of either rate: it
renormalises the forward constant by a factor independent of the loading. So
when the salt side is swept (c_L2 = c_LF = 0, hence a0 = b0 = 0) the atomic
fluxes are exactly

    J_R = 2 w_rec = 2 A c**2,    J_F = w_F = B1 c

and the apparent exponent of Eq. (apparent_exponent) read off the total flux is
exactly Eq. (n_of_B),

    n = dln J / dln c = (2 J_R + J_F) / (J_R + J_F) = (2 + B) / (1 + B)

with B = J_F/J_R the branching ratio of Eq. (branching). The reverse terms are
active throughout; they simply do not move the exponent.

Second, the same exponent read off the salt-side *inventory* rather than the
flux is not the same number. The total interfacial salt-side hydrogen is
2 c_H2|G + c_HF|G = J_R R_2 + J_F R_F, so its logarithmic slope is
(2 + B R_F/R_2)/(1 + B R_F/R_2): B is replaced by B R_HF/R_H2, and the two
readouts coincide only when the two carriers are transported alike. This is the
last paragraph of the analytical appendix.
"""

import numpy as np


def _positive_root(a, b, c):
    """Root of a x**2 + b x + c = 0 that is positive and stays finite as a -> 0.

    Here b > 0 and c < 0 for every case, so the discriminant exceeds b**2 and
    exactly one root is positive. The citardauq form below evaluates it without
    the cancellation that (-b + sqrt(disc)) suffers when a is small, and tends
    to -c/b as a -> 0.
    """
    disc = b**2 - 4 * a * c
    if disc < 0:
        raise ValueError(f"no real steady state: discriminant {disc:.3e} < 0")
    return 2 * c / (-b - np.sqrt(disc))


def solution(c_0, c_L2, c_LF, D_m, D_H2, D_HF, L_m, L_s, kr_plus, kr_minus, kf_plus, kf_minus, a_F):
    """Steady state of the two-channel interface.

    Args:
        c_0: atomic concentration imposed on the outer metal face
        c_L2, c_LF: H2 and HF concentrations imposed on the outer liquid face
        D_m: diffusivity of atomic H in the metal
        D_H2, D_HF: diffusivities of the two carriers in the liquid
        L_m, L_s: slab thicknesses
        kr_plus, kr_minus: constants of the recombination channel
        kf_plus, kf_minus: constants of the fluorination channel
        a_F: effective fluoride activity, folded into the forward constant

    Returns:
        (c_m, c_H2, c_HF, w_rec, w_F) at the interface. The atomic flux leaving
        the metal is 2 w_rec + w_F.
    """
    R_m = L_m / D_m
    R_2 = L_s / D_H2
    R_F = L_s / D_HF

    # each channel rate as an explicit function of the metal-side trace
    A = kr_plus / (1 + kr_minus * R_2)
    a0 = kr_minus * c_L2 / (1 + kr_minus * R_2)
    B1 = kf_plus * a_F / (1 + kf_minus * R_F)
    b0 = kf_minus * c_LF / (1 + kf_minus * R_F)

    # metal balance c = c_0 - (2 w_rec + w_F) R_m
    c_m = _positive_root(
        2 * A * R_m,
        1 + B1 * R_m,
        -(c_0 + R_m * (2 * a0 + b0)),
    )

    w_rec = A * c_m**2 - a0
    w_F = B1 * c_m - b0

    return c_m, c_L2 + w_rec * R_2, c_LF + w_F * R_F, w_rec, w_F


def branching(w_rec, w_F):
    """Branching ratio of Eq. (branching), B = w_F / (2 w_rec): the atomic flux
    carried by fluorination over that carried by recombination."""
    return w_F / (2 * w_rec)


def apparent_exponent(B):
    """Eq. (n_of_B): n = (2 + B)/(1 + B), the logarithmic slope of the total
    atomic flux against the interfacial loading. Falls monotonically from 2 at
    B -> 0 (recombination only, Sieverts/Henry) to 1 at B -> infinity
    (fluorination only), through 3/2 at B = 1."""
    return (2 + B) / (1 + B)
