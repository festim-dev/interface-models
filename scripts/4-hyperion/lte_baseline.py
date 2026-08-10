"""LTE reference transient for the Ni/FLiBe operating point of parameters.py.

This is the baseline the kinetic interface models have to reproduce in the
fast-kinetics limit, and it is the curve every later figure in Sec. 5 is read
against. The interface is treated the way a macroscopic hydrogen transport code
treats it: an algebraic per-species constraint, here Sieverts on the metal
against Henry in the salt, imposed by `F.Interface` with the penalty method.
That is the same condition `para_1d.py` in the HYPERION modelling repo uses, made
transient.

A single species carries hydrogen on both sides, so the salt-side variable is
an atomic concentration and the Henry constant is the atomic one. Model 2 gives
the salt a molecular carrier instead and needs the molecular constant; the two
differ by the stoichiometric factor and parameters.py holds both. Keeping the
conventions straight is what makes the fast-kinetics limit land on this curve
rather than a factor of two away from it.

Verification: at the final time the computed interfacial concentrations and the
downstream flux are compared with the analytical solution in
`parameters.lte_steady_state`. The steady profiles are piecewise linear and so
nodally exact on P1, which leaves the interface condition alone under test.

Produces:

  - lte_baseline.csv : time, downstream flux, interfacial concentrations
  - lte_baseline.pdf : the transient and the steady profile
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np

import festim as F

import parameters as p

# Penalty weight. The constraint is (u_m/K_S)^2 - u_s/K_H, whose terms are of
# order 1e4 here, so the penalty has to be large enough to dominate the
# diffusive residual. `para_1d.py` uses 1e20 for the steady problem and the
# verification block below checks it is doing its job rather than assuming it.
PENALTY_TERM = 1e20

FINAL_TIME = 1.0e5  # s, about five salt diffusion times
INITIAL_STEP = 10.0  # s
# The time integration is backward Euler and so first order, and the time lag
# measured in sweep_damkohler.py is sensitive to the early transient: the
# absolute lag is still drifting at 2000 s steps. 250 s costs little here and
# the ratio of two lags, which is what the sweep reports, is converged to a few
# parts in ten thousand by then.
MAX_STEP = 250.0  # s


def graded_vertices(n_metal=160, n_salt=220, bias=2.5):
    """1D mesh refined towards the interface from both sides.

    The transient front sits in the salt for most of the run, but the interface
    is where the constraint acts and where the kinetic models will put a
    boundary layer, so both sides are graded towards x = L_NI.
    """
    s = np.linspace(0.0, 1.0, n_metal)
    metal = p.L_NI * (1.0 - (1.0 - s) ** bias)

    s = np.linspace(0.0, 1.0, n_salt)
    salt = p.L_NI + p.L_SALT * (s**bias)

    return np.unique(np.concatenate([metal, salt]))


def build_model(temperature=p.TEMPERATURE, p_up=p.P_UP, p_down=p.P_DOWN):
    """Two-slab transient permeation problem with LTE at the interface."""
    metal = F.Material(
        D_0=p.D_0_NI,
        E_D=p.E_D_NI,
        K_S_0=p.K_S_0_NI,
        E_K_S=p.E_K_S_NI,
        solubility_law="sievert",
    )
    salt = F.Material(
        D_0=p.D_0_SALT,
        E_D=p.E_D_SALT,
        K_S_0=p.K_H_0_SALT,
        E_K_S=p.E_K_H_SALT,
        solubility_law="henry",
    )

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

    H = F.Species("H", subdomains=[vol_metal, vol_salt])
    model.species = [H]
    model.temperature = temperature

    model.interfaces = [
        F.Interface(
            id=5,
            subdomains=[vol_metal, vol_salt],
            penalty_term=PENALTY_TERM,
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
            species=H,
            value=p.K_H_atomic(temperature) * p_down,
        ),
    ]

    flux = F.SurfaceFlux(field=H, surface=downstream)
    profile_metal = F.Profile1DExport(field=H, subdomain=vol_metal)
    profile_salt = F.Profile1DExport(field=H, subdomain=vol_salt)
    model.exports = [flux, profile_metal, profile_salt]

    # The residual is dominated by the penalty term, which is larger than the
    # diffusive terms by construction, so a tolerance chosen from the size of
    # the concentrations is far too loose: atol=1e15 lets Newton stop early and
    # the final state misses the analytical steady state by a few parts per
    # thousand, which looks like a discretisation error and is not one. These
    # are the values para_1d.py uses. Tightening further (atol=1e8) stops SNES
    # converging. This is a concrete instance of the interface-residual scaling
    # listed as open in Sec. "Numerical considerations".
    model.settings = F.Settings(
        atol=1e10,
        rtol=1e-12,
        max_iterations=100,
        transient=True,
        final_time=FINAL_TIME,
        stepsize=F.Stepsize(
            initial_value=INITIAL_STEP,
            growth_factor=1.2,
            cutback_factor=0.8,
            target_nb_iterations=4,
            max_stepsize=MAX_STEP,
        ),
    )

    return model, flux, profile_metal, profile_salt


def run(**kwargs):
    """Solve and return (times, flux, c_metal_gamma, c_salt_gamma, model)."""
    model, flux, profile_metal, profile_salt = build_model(**kwargs)
    model.initialise()
    model.run()

    times = np.array(profile_metal.t)
    # SurfaceFlux gives the outward flux; permeation out of the salt is
    # negative on that convention, so take the magnitude
    flux_history = np.abs(np.array(flux.data))
    c_metal = np.array([row[-1] for row in profile_metal.data])
    c_salt = np.array([row[0] for row in profile_salt.data])

    return times, flux_history, c_metal, c_salt, model


def verify(times, flux_history, c_metal, c_salt):
    """Compare the final state with the analytical LTE steady state."""
    c_0, c_gamma_ref, c_salt_ref, flux_ref = p.lte_steady_state()

    errors = {
        "c_metal|Gamma": abs(c_metal[-1] / c_gamma_ref - 1.0),
        "c_salt|Gamma": abs(c_salt[-1] / c_salt_ref - 1.0),
        "flux": abs(flux_history[-1] / flux_ref - 1.0),
    }

    print(f"\nfinal time                 {times[-1]:.4e} s")
    print(f"c_metal|Gamma  computed    {c_metal[-1]:.6e} m^-3")
    print(f"               analytical solution {c_gamma_ref:.6e} m^-3")
    print(f"c_salt|Gamma   computed    {c_salt[-1]:.6e} m^-3")
    print(f"               analytical solution {c_salt_ref:.6e} m^-3")
    print(f"flux           computed    {flux_history[-1]:.6e} m^-2 s^-1")
    print(f"               analytical solution {flux_ref:.6e} m^-2 s^-1")
    print()
    for name, err in errors.items():
        print(f"relative error {name:<16s} {err:.3e}")

    # the interface constraint itself, which is what the penalty term enforces
    lhs = (c_metal[-1] / p.K_S_nickel()) ** 2
    rhs = c_salt[-1] / p.K_H_atomic()
    print(f"\ninterface constraint (u_m/K_S)^2 = {lhs:.6e} Pa")
    print(f"                     u_s/K_H     = {rhs:.6e} Pa")
    print(f"                     mismatch    = {abs(lhs / rhs - 1.0):.3e}")

    return errors


def plot(times, flux_history, model, filename):
    """(a) the permeation transient, (b) the steady profile across both slabs."""
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    _, _, _, flux_ref = p.lte_steady_state()
    tau_metal, tau_salt = p.diffusion_times()
    GUIDE = "0.45"

    fig, (ax_transient, ax_profile) = plt.subplots(1, 2, figsize=(9.0, 3.8))

    ax_transient.plot(times, flux_history, color="C0", linewidth=1.6)
    ax_transient.axhline(flux_ref, color=GUIDE, linestyle="--", linewidth=1)
    ax_transient.annotate(
        "analytical steady state",
        xy=(times[-1], flux_ref),
        xytext=(-6, -14),
        textcoords="offset points",
        ha="right",
        fontsize=8,
        color=GUIDE,
        weight="bold",
    )
    for tau, label in ((tau_metal, r"$\tau_\mathrm{Ni}$"), (tau_salt, r"$\tau_\mathrm{salt}$")):
        ax_transient.axvline(tau, color="C2", linestyle=":", linewidth=1)
        ax_transient.annotate(
            label,
            xy=(tau, flux_ref * 0.08),
            xytext=(4, 0),
            textcoords="offset points",
            fontsize=9,
            color="C2",
            weight="bold",
        )
    ax_transient.set_xlabel("time [s]")
    ax_transient.set_ylabel(r"downstream flux [m$^{-2}$ s$^{-1}$]")
    ax_transient.set_xlim(0, times[-1])
    ax_transient.set_ylim(bottom=0)
    ax_transient.set_title("(a) permeation transient", loc="left", fontsize=10)

    profile_metal, profile_salt = model.exports[1], model.exports[2]
    for export, colour, label in (
        (profile_metal, "C0", "Ni"),
        (profile_salt, "C1", "FLiBe"),
    ):
        ax_profile.plot(
            np.array(export.x) * 1e3, export.data[-1], color=colour, linewidth=1.6
        )
        ax_profile.annotate(
            label,
            xy=(np.mean(export.x) * 1e3, np.mean(export.data[-1])),
            fontsize=9,
            color=colour,
            weight="bold",
        )

    ax_profile.axvline(p.L_NI * 1e3, color="#f46036", linewidth=1.5)
    ax_profile.annotate(
        r"$\Gamma$",
        xy=(p.L_NI * 1e3, 0),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=11,
        color="#f46036",
        weight="bold",
    )
    ax_profile.set_xlabel("$x$ [mm]")
    ax_profile.set_ylabel(r"$c$ [m$^{-3}$]")
    ax_profile.set_ylim(bottom=0)
    ax_profile.set_title("(b) steady profile", loc="left", fontsize=10)

    fig.tight_layout()
    fig.savefig(filename)
    print(f"\nwrote {filename}")


if __name__ == "__main__":
    times, flux_history, c_metal, c_salt, model = run()
    verify(times, flux_history, c_metal, c_salt)

    np.savetxt(
        "lte_baseline.csv",
        np.column_stack([times, flux_history, c_metal, c_salt]),
        delimiter=",",
        header="t,flux,c_metal_gamma,c_salt_gamma",
        comments="",
    )
    plot(times, flux_history, model, "lte_baseline.pdf")
