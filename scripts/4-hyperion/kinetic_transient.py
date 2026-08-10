"""Model 2 on the Ni/FLiBe operating point: recombination into a molecular carrier.

The same problem as lte_baseline.py, with the algebraic interface constraint
replaced by one reversible mass-action channel,

    R:   2 H(m)  <->  H2(s) ,   w = k_+ c_H|G^2 - k_- c_H2|G ,

declared as a single `InterfaceReaction` with `reactants=[H] * 2`. The
repetition puts the square on the rate and the stoichiometric factor two on the
metal-side flux, so the atomic flux leaving the metal is 2w while the salt gains
one molecule per event.

Detailed balance fixes k_+/k_- = K_H_molecular / K_S^2 from the same
thermodynamics LTE uses, which is what makes this a strict generalisation and
not a competing model. FESTIM does not enforce it, so parameters.py imposes it
and `model2_steady_state` below reduces exactly to the LTE analytical solution in the
fast-kinetics limit. That reduction is analytic, not numerical: the coefficient
A of the quadratic goes over to the LTE one term by term as k_+ grows at fixed
ratio. The numbers printed by this script check that FESTIM agrees.

Note the two Henry conventions. The salt-side unknown here counts H2 molecules,
so the constant that appears is K_H_molecular, half the atomic one the LTE
baseline uses. Getting that wrong moves the fast-kinetics limit by a factor of
two and is the easiest way to conclude, wrongly, that the kinetic model does not
recover LTE.

Produces, when run directly:

  - a table of FESTIM against the analytical solution across six decades of Da
  - kinetic_transient.csv
"""

import matplotlib

matplotlib.use("Agg")

import numpy as np

import festim as F

import parameters as p
from lte_baseline import FINAL_TIME, INITIAL_STEP, MAX_STEP, graded_vertices


def model2_steady_state(k_plus, k_minus=None, temperature=p.TEMPERATURE, p_up=p.P_UP):
    """Analytical steady state of the single recombination channel.

    Both bulk profiles are linear at steady state, so eliminating them leaves
    an algebraic problem in the interfacial loading. The salt carries molecules
    at rate w, so c_H2|G = w L_s / D_s with the outer face swept, and

        w = k_+ c_m^2 / (1 + k_- L_s / D_s) ,

    while the metal delivers 2w = D_m (c_0 - c_m) / L_m. Together,

        A c_m^2 + c_m - c_0 = 0 ,
        A = 2 k_+ L_m / [ D_m (1 + k_- L_s / D_s) ] ,

    and the positive root is taken. As k_+ grows at fixed k_+/k_-, A tends to
    (K_H_atomic / K_S^2) L_m D_s / (L_s D_m), which is exactly the coefficient
    of the LTE quadratic in parameters.lte_steady_state. The two analytical solutions
    are therefore the same equation in that limit, term by term.

    Returns (c_0, c_metal|G, c_H2|G, J_atomic).
    """
    if k_minus is None:
        k_minus = p.detailed_balance_k_minus(k_plus, temperature)

    c_0 = p.upstream_concentration(temperature, p_up)
    D_m, D_s = p.D_nickel(temperature), p.D_salt(temperature)

    A = 2 * k_plus * p.L_NI / (D_m * (1 + k_minus * p.L_SALT / D_s))
    c_metal = (-1 + np.sqrt(1 + 4 * A * c_0)) / (2 * A)

    w = k_plus * c_metal**2 / (1 + k_minus * p.L_SALT / D_s)
    c_h2 = w * p.L_SALT / D_s

    return c_0, c_metal, c_h2, 2 * w


def build_model(
    k_plus,
    k_minus=None,
    temperature=p.TEMPERATURE,
    p_up=p.P_UP,
    p_down=p.P_DOWN,
    transient=True,
    final_time=FINAL_TIME,
):
    """Two-slab problem with one reversible channel on the interface.

    The metal carries atomic H and the salt carries molecular H2, so unlike the
    LTE baseline there is no single species spanning both subdomains and no
    constraint tying the two traces together. Everything the interface does is
    in the reaction term.
    """
    if k_minus is None:
        k_minus = p.detailed_balance_k_minus(k_plus, temperature)

    metal = F.Material(D_0=p.D_0_NI, E_D=p.E_D_NI)
    salt = F.Material(D_0=p.D_0_SALT, E_D=p.E_D_SALT)

    x_gamma = p.L_NI
    x_end = p.L_NI + p.L_SALT

    vol_metal = F.VolumeSubdomain1D(id=1, borders=[0.0, x_gamma], material=metal)
    vol_salt = F.VolumeSubdomain1D(id=2, borders=[x_gamma, x_end], material=salt)
    upstream = F.SurfaceSubdomain1D(id=3, x=0.0)
    downstream = F.SurfaceSubdomain1D(id=4, x=x_end)

    model = F.HydrogenTransportProblemDiscontinuous()
    model.subdomains = [vol_metal, vol_salt, upstream, downstream]
    model.surface_to_volume = {upstream: vol_metal, downstream: vol_salt}
    model.mesh = F.Mesh1D(graded_vertices())

    H = F.Species("H", subdomains=[vol_metal])
    H2 = F.Species("H2", subdomains=[vol_salt])
    model.species = [H, H2]
    model.temperature = temperature

    model.interfaces = [
        F.InterfaceReaction(
            id=5,
            subdomains=[vol_metal, vol_salt],
            k_plus=k_plus,
            k_minus=k_minus,
            reactants=[H] * 2,
            products=[H2],
        )
    ]

    model.boundary_conditions = [
        F.SievertsBC(
            subdomain=upstream,
            species=H,
            pressure=p_up,
            S_0=p.K_S_0_NI,
            E_S=p.E_K_S_NI,
        ),
        F.FixedConcentrationBC(
            subdomain=downstream,
            species=H2,
            value=p.K_H_molecular(temperature) * p_down,
        ),
    ]

    flux = F.SurfaceFlux(field=H2, surface=downstream)
    profile_metal = F.Profile1DExport(field=H, subdomain=vol_metal)
    profile_salt = F.Profile1DExport(field=H2, subdomain=vol_salt)
    model.exports = [flux, profile_metal, profile_salt]

    # Same tolerances as the LTE baseline and for the same reason: the residual
    # scale is set by the interface term, not by the concentrations.
    model.settings = F.Settings(
        atol=1e10,
        rtol=1e-12,
        max_iterations=100,
        transient=transient,
        final_time=final_time if transient else None,
        stepsize=F.Stepsize(
            initial_value=INITIAL_STEP,
            growth_factor=1.2,
            cutback_factor=0.8,
            target_nb_iterations=4,
            max_stepsize=MAX_STEP,
        )
        if transient
        else None,
    )

    return model, flux, profile_metal, profile_salt


def run(k_plus, **kwargs):
    """Solve and return (times, atomic flux, c_metal|G, c_H2|G, model).

    The exported surface flux counts molecules, so it is doubled to give the
    atomic flux the LTE baseline reports and a permeation experiment measures.
    """
    model, flux, profile_metal, profile_salt = build_model(k_plus, **kwargs)
    model.initialise()
    model.run()

    times = np.array(profile_metal.t)
    flux_history = 2.0 * np.abs(np.array(flux.data))
    c_metal = np.array([row[-1] for row in profile_metal.data])
    c_h2 = np.array([row[0] for row in profile_salt.data])

    return times, flux_history, c_metal, c_h2, model


if __name__ == "__main__":
    _, c_gamma_lte, _, flux_lte = p.lte_steady_state()

    print("Detailed balance check: one channel, k_-/k_+ from thermodynamics.")
    print("As Da grows the kinetic steady state must approach the LTE one.\n")
    print(f"LTE analytical solution: c_m|G = {c_gamma_lte:.6e} m^-3, "
          f"J = {flux_lte:.6e} m^-2 s^-1\n")
    print(f"{'Da':>10s} {'c_m|G':>14s} {'J':>14s} {'J/J_LTE':>10s} "
          f"{'err vs analytical':>20s}")

    rows = []
    for damkohler in np.logspace(-2, 4, 7):
        k_plus = p.k_plus_from_damkohler(damkohler)
        times, flux_history, c_metal, c_h2, _ = run(k_plus)

        _, c_ref, _, flux_ref = model2_steady_state(k_plus)
        err = max(
            abs(c_metal[-1] / c_ref - 1.0), abs(flux_history[-1] / flux_ref - 1.0)
        )

        print(f"{damkohler:10.2e} {c_metal[-1]:14.6e} {flux_history[-1]:14.6e} "
              f"{flux_history[-1] / flux_lte:10.6f} {err:20.3e}")

        rows.append(
            [damkohler, k_plus, p.detailed_balance_k_minus(k_plus),
             c_metal[-1], c_h2[-1], flux_history[-1],
             c_ref, flux_ref, flux_history[-1] / flux_lte, err]
        )

    np.savetxt(
        "kinetic_transient.csv",
        np.array(rows),
        delimiter=",",
        header="Da,k_plus,k_minus,c_metal,c_H2,J,c_metal_ref,J_ref,J_over_J_lte,err",
        comments="",
    )
    print("\nwrote kinetic_transient.csv")
