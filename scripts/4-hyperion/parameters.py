"""Operating point for Sec. "Application: the HYPERION experiment".

Every dimensional number used in Sec. 5 is defined here and nowhere else, so
the values quoted in the paper have one place to be checked against. Nothing in
this module is fitted and nothing is compared with measurement: HYPERION
supplies a realistic geometry and operating window, and the simulations that
import this file are a demonstration of what the kinetic interface framework
predicts there.

Provenance
----------
Geometry and operating conditions are read off the HYPERION modelling repo,
github.com/festim-dev/hyperion at commit 8f61a8c:

  L_NI, RADIUS          para_1d.py, `run_all_cases_1d` defaults
  L_SALT                para_1d.py, `L_FLIBE_BY_TEMP_C`
  vessel coordinates    mesh.py, `generate_mesh`
  P_UP                  exp_data.py, `swap_conditions`, a measured upstream
                        pressure at 500 C

Transport properties are Arrhenius fits from h_transport_materials, hard-coded
here rather than imported so that the paper repo pins the numbers and does not
move when the database is updated. One source per material:

  nickel    Louthan et al. 1975
  flibe     Calderoni et al. 2008   (`calderoni_measurement_2008` in the bib)

Isotope mass effects on transport are neglected: D and K are taken identical
for H and T, so entries measured with different isotopes are used
interchangeably. This is consistent with the mass-independent limit used for
the isotopologue rate constants in Sec. "Multiple isotopes", where
k_HT+ = 2 k_HH+ is a degeneracy factor and not a mass effect, and it means any
isotope separation the model produces is unambiguously interfacial.

Units
-----
SI throughout, concentrations in particles per cubic metre, pressures in Pa,
energies in eV.

The salt-side Henry constant needs a convention, and which one is meant is
exactly the question of Sec. "A solubility constant has no law-independent
units". K_H_ATOMIC counts hydrogen atoms; K_H_MOLECULAR counts H2 molecules and
is half of it. The LTE baseline carries a single species and uses the atomic
form; Model 2 carries a molecular carrier and uses the molecular form. Reported
values do not say which is meant, and the factor of two between them is the
smallest instance of the ambiguity the section is about.
"""

import numpy as np

# ── physical constants ───────────────────────────────────────────────────────

K_B = 8.617333262e-5  # Boltzmann constant [eV/K]

# ── geometry ─────────────────────────────────────────────────────────────────

L_NI = 2.032e-3  # Ni membrane thickness [m] (0.080 in)
L_SALT = 5.139858e-3  # FLiBe layer thickness at 500 C [m]
RADIUS = 3.07 * 0.0254 / 2  # membrane radius [m], 3.07 in diameter

MEMBRANE_AREA = np.pi * RADIUS**2  # [m^2]
SIDEWALL_AREA = 2 * np.pi * RADIUS * L_SALT  # wetted sidewall [m^2]

# Vessel coordinates for the sketch, from mesh.py. Radii are (inner, outer) and
# the layers run bottom to top along the axis. The 1D model keeps the membrane
# and the pool and drops the rest; the sidewall it drops carries a permeation
# path in parallel with the modelled one, which is why SIDEWALL_AREA is
# reported alongside MEMBRANE_AREA.
R_INNER = 0.039  # [m]
R_OUTER = 0.041  # [m]
VESSEL_LAYERS = (  # (name, y_bottom [m], y_top [m], kind)
    ("bottom cap", 0.000, 0.002, "nickel"),
    ("upstream plenum", 0.002, 0.022, "gas"),
    ("membrane", 0.022, 0.024, "nickel"),
    ("FLiBe pool", 0.024, 0.02914, "salt"),
    ("cover gas", 0.02914, 0.1091, "gas"),
    ("top cap", 0.1091, 0.1111, "nickel"),
)

# mesh.py rounds the membrane to 2.0 mm while para_1d.py uses 2.032 mm. The 1D
# model uses L_NI; the sketch is drawn to the mesh.py coordinates above.

# ── operating point ──────────────────────────────────────────────────────────

TEMPERATURE = 773.15  # [K], 500 C, the low end of the HYPERION window
P_UP = 1.31e5  # upstream H2 pressure [Pa]
P_DOWN = 0.0  # swept free surface, idealised

# The rate constants carry no Arrhenius dependence on the interface-flux branch
# of FESTIM, so Sec. 5 is run at a single temperature and says so.

# ── transport properties ─────────────────────────────────────────────────────

D_0_NI = 7.0e-7  # [m^2/s]      Louthan 1975
E_D_NI = 0.40938865142235586  # [eV]
K_S_0_NI = 3.312177418e23  # [m^-3 Pa^-1/2]  Louthan 1975
E_K_S_NI = 0.16375546056894236  # [eV]

D_0_SALT = 9.3e-7  # [m^2/s]      Calderoni 2008
E_D_SALT = 0.4352993255630113  # [eV]
K_H_0_SALT = 4.7574912004e22  # [m^-3 Pa^-1]   Calderoni 2008
E_K_H_SALT = 0.3627494379691761  # [eV]


def arrhenius(pre_exponential, activation_energy, temperature=TEMPERATURE):
    """Arrhenius law, energies in eV."""
    return pre_exponential * np.exp(-activation_energy / (K_B * temperature))


def D_nickel(temperature=TEMPERATURE):
    """Diffusivity of hydrogen in nickel [m^2/s]."""
    return arrhenius(D_0_NI, E_D_NI, temperature)


def K_S_nickel(temperature=TEMPERATURE):
    """Sieverts constant of hydrogen in nickel [m^-3 Pa^-1/2]."""
    return arrhenius(K_S_0_NI, E_K_S_NI, temperature)


def D_salt(temperature=TEMPERATURE):
    """Diffusivity of the hydrogen carrier in FLiBe [m^2/s]."""
    return arrhenius(D_0_SALT, E_D_SALT, temperature)


def K_H_atomic(temperature=TEMPERATURE):
    """Henry constant of FLiBe counting hydrogen atoms [m^-3 Pa^-1]."""
    return arrhenius(K_H_0_SALT, E_K_H_SALT, temperature)


def K_H_molecular(temperature=TEMPERATURE):
    """Henry constant of FLiBe counting H2 molecules [m^-3 Pa^-1].

    Half the atomic form. See the units note in the module docstring: the two
    differ by the stoichiometric factor and reported values do not say which
    convention they follow.
    """
    return 0.5 * K_H_atomic(temperature)


# ── derived scales ───────────────────────────────────────────────────────────


def diffusion_times(temperature=TEMPERATURE):
    """(tau_metal, tau_salt) in seconds, the two bulk diffusion times."""
    return L_NI**2 / D_nickel(temperature), L_SALT**2 / D_salt(temperature)


def upstream_concentration(temperature=TEMPERATURE, pressure=P_UP):
    """Sieverts loading on the upstream face of the membrane [m^-3]."""
    return K_S_nickel(temperature) * np.sqrt(pressure)


def lte_steady_state(temperature=TEMPERATURE, p_up=P_UP):
    """Analytical LTE steady state of the two-slab problem, atomic convention.

    Both profiles are linear at steady state. Continuity of the atomic flux and
    equality of the equilibrium pressure across the interface give

        D_m (c_0 - c_G) / L_m = D_s K_H (c_G / K_S)^2 / L_s ,

    with the downstream face swept, and the positive root is taken. Returns
    (c_upstream, c_interface_metal, c_interface_salt, J_atomic), the flux per
    unit area in particles per square metre per second.

    This is the number lte_baseline.py has to reproduce, and through detailed
    balance it is also the large-Damkoehler limit the kinetic models have to
    reproduce.
    """
    c_0 = upstream_concentration(temperature, p_up)
    K_S = K_S_nickel(temperature)

    a = D_salt(temperature) * K_H_atomic(temperature) / (L_SALT * K_S**2)
    b = D_nickel(temperature) / L_NI

    c_gamma = (-b + np.sqrt(b**2 + 4 * a * b * c_0)) / (2 * a)
    c_salt = K_H_atomic(temperature) * (c_gamma / K_S) ** 2

    return c_0, c_gamma, c_salt, b * (c_0 - c_gamma)


def detailed_balance_k_minus(k_plus, temperature=TEMPERATURE):
    """Reverse constant of the recombination channel 2 H(m) <-> H2(s).

    Detailed balance ties the ratio to the same thermodynamics LTE uses,

        k_+ / k_- = K_H_molecular / K_S^2 ,

    which is what makes the kinetic framework a strict generalisation and not a
    competing model. FESTIM does not enforce this, so it is imposed here.
    """
    ratio = K_H_molecular(temperature) / K_S_nickel(temperature) ** 2
    return k_plus / ratio


def damkohler_recombination(k_plus, c_ref=None, temperature=TEMPERATURE, p_up=P_UP):
    """Damkoehler number of the recombination channel.

    The channel is second order, so `k_plus` is a m^4/s and not a velocity, and
    it has to be made dimensionless with a loading as well as with the
    metal-side transport scale. Following Eq. (damkohler) of the paper the
    channel is linearised about the interfacial state, which gives the exchange
    velocity 2 k_+ c_m|G, and the reference loading is the upstream Sieverts
    value c* = K_S sqrt(P_up):

        Da = 2 k_+ c* L_m / D_m .

    The same convention runs through the paper, in
    2-higher-order-reactions/verification_model2_lte_limit.py as well as here:
    c* is known from the boundary condition before anything is solved, so Da is
    an input to a sweep, and the value actually attained at the interface is
    reported as a diagnostic.
    """
    if c_ref is None:
        c_ref = upstream_concentration(temperature, p_up)
    return 2 * k_plus * c_ref * L_NI / D_nickel(temperature)


def k_plus_from_damkohler(damkohler, temperature=TEMPERATURE, p_up=P_UP):
    """Forward recombination constant giving a target Damkoehler number."""
    c_ref = upstream_concentration(temperature, p_up)
    return damkohler * D_nickel(temperature) / (2 * c_ref * L_NI)


if __name__ == "__main__":
    tau_m, tau_s = diffusion_times()
    c_0, c_gamma, c_salt, flux = lte_steady_state()

    print(f"temperature            {TEMPERATURE:.2f} K ({TEMPERATURE - 273.15:.0f} C)")
    print(f"upstream pressure      {P_UP:.3e} Pa")
    print()
    print(f"D_nickel               {D_nickel():.4e} m^2/s")
    print(f"K_S_nickel             {K_S_nickel():.4e} m^-3 Pa^-1/2")
    print(f"D_salt                 {D_salt():.4e} m^2/s")
    print(f"K_H_atomic             {K_H_atomic():.4e} m^-3 Pa^-1")
    print(f"K_H_molecular          {K_H_molecular():.4e} m^-3 Pa^-1")
    print()
    print(f"tau_metal              {tau_m:.4e} s")
    print(f"tau_salt               {tau_s:.4e} s   ({tau_s / tau_m:.1f} x tau_metal)")
    print()
    print(f"c upstream             {c_0:.4e} m^-3")
    print(f"c metal at interface   {c_gamma:.4e} m^-3  ({100 * c_gamma / c_0:.1f} % of upstream)")
    print(f"c salt  at interface   {c_salt:.4e} m^-3")
    print(f"LTE steady flux        {flux:.4e} m^-2 s^-1")
    print(f"                       {flux * MEMBRANE_AREA:.4e} s^-1 over {MEMBRANE_AREA:.4e} m^2")
    print()
    print(f"membrane area          {MEMBRANE_AREA:.4e} m^2")
    print(f"wetted sidewall area   {SIDEWALL_AREA:.4e} m^2  "
          f"({100 * SIDEWALL_AREA / (SIDEWALL_AREA + MEMBRANE_AREA):.0f} % of wetted Ni)")
    print()
    print(f"k_+/k_- (detailed balance)  {K_H_molecular() / K_S_nickel() ** 2:.4e} m^3")
