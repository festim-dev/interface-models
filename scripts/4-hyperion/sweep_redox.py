"""Experiment B: two channels on the interface, sweeping the branching ratio.

Model 3. The recombination channel of kinetic_transient.py is joined by a
fluorination channel on the same interface,

    R:   2 H(m)  <->  H2(s)
    F:     H(m)  <->  HF(s)

declared as two `InterfaceReaction` objects sharing an interface id, so their
residual contributions add. The bulk is untouched: same diffusivities, same
solubilities, same geometry, same upstream pressure as every other run in this
directory.

Parameterisation, and why it is not a solubility
------------------------------------------------
The fluorination channel has no measured thermodynamics to anchor to. There is
no Henry constant for HF in FLiBe in h_transport_materials and none in the
paper's bibliography, and by the argument of Sec. "A solubility constant has no
law-independent units" there could not be one that is independent of an assumed
carrier. So the channel is not given a solubility. It is driven by the two
dimensionless groups the paper already defines:

    Da   = 2 k_r+ c* L_m / D_m           the recombination Damkoehler number
    B    = w_F / (2 w_rec)               the branching ratio, Eq. (branching)

both built at the same reference loading c* = K_S sqrt(P_up), the upstream
Sieverts value, as in parameters.damkohler_recombination. Far from equilibrium
B -> k_f+ a_F / (2 k_r+ c_m|G), so a target B at that reference fixes the
forward constant of the F channel,

    k_f+ a_F = 2 B k_r+ c* ,

and the Damkoehler number of the F channel is then not free: Da_F = B Da.
a_F never appears on its own, only folded into the forward constant, which is
also all the implementation supports.

The reverse constant needs one more choice. In the steady state the reverse
term enters only through the dimensionless combination rho = k_- L_s / D, since
the swept downstream face gives c|G = w L_s / D for each carrier, so

    w = k_+ (reactants) / (1 + rho) .

rho is therefore the natural measure of how much the reverse term renormalises a
channel, and detailed balance fixes it for R. For F we set rho_F = rho_R, so
both channels sit at the same distance from their own equilibrium, and we state
it as the assumption it is. Because rho only renormalises the forward constant,
this choice moves the *nominal* B and not the physics; the branching ratio is
therefore measured from the computed rates afterwards, as in
scripts/3-competing-channels/verification_model3_exponent.py, and everything is
plotted against the measured value.

Both carriers are given the same diffusivity, for want of a measurement for HF.
That is worth remembering when reading the flux split: it means the split
reported here is set at the interface and not by differential transport away
from it.

Produces:

  - sweep_redox.csv
  - sweep_redox.pdf
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

import festim as F

import parameters as p
import lte_baseline
from lte_baseline import FINAL_TIME, INITIAL_STEP, MAX_STEP, graded_vertices
from sweep_damkohler import time_lag

DAMKOHLER = 10.0  # recombination channel, held fixed across the sweep
ALL_B = np.logspace(-2, 2, 9)  # target branching ratio

GUIDE = "0.45"


def channel_constants(damkohler=DAMKOHLER, branching=1.0, temperature=p.TEMPERATURE):
    """Forward and reverse constants of both channels for a target (Da, B).

    Returns (k_r_plus, k_r_minus, k_f_plus, k_f_minus). k_r_minus comes from
    detailed balance; k_f_minus from rho_F = rho_R as set out in the module
    docstring.
    """
    c_ref = p.upstream_concentration(temperature)

    k_r_plus = p.k_plus_from_damkohler(damkohler, temperature)
    k_r_minus = p.detailed_balance_k_minus(k_r_plus, temperature)

    k_f_plus = 2.0 * branching * k_r_plus * c_ref

    # rho = k_minus L_s / D, equal for the two channels; both carriers share D
    k_f_minus = k_r_minus

    return k_r_plus, k_r_minus, k_f_plus, k_f_minus


def model3_steady_state(damkohler=DAMKOHLER, branching=1.0, temperature=p.TEMPERATURE,
                        p_up=p.P_UP):
    """Closed-form steady state of the two-channel interface.

    Both bulk profiles are linear at steady state and the downstream face is
    swept, so each carrier satisfies c|G = w L_s / D and each channel collapses
    to its forward term with the reverse folded in,

        w_rec = k_r+ c_m^2 / (1 + k_r- L_s / D_s) = A c_m^2 ,
        w_F   = k_f+ c_m   / (1 + k_f- L_s / D_s) = C c_m .

    The metal delivers 2 w_rec + w_F = D_m (c_0 - c_m) / L_m, so the
    interfacial loading is the positive root of

        2 A c_m^2 + (C + D_m/L_m) c_m - D_m c_0 / L_m = 0 ,

    and the atomic flux follows. Setting C = 0 recovers the Model 2 form in
    kinetic_transient.py, which in turn recovers the LTE one; the three closed
    forms are nested.

    Returns (c_metal|G, J_atomic, HF fraction of the atomic flux).
    """
    k_r_plus, k_r_minus, k_f_plus, k_f_minus = channel_constants(
        damkohler, branching, temperature
    )
    c_0 = p.upstream_concentration(temperature, p_up)
    D_m, D_s = p.D_nickel(temperature), p.D_salt(temperature)

    A = k_r_plus / (1.0 + k_r_minus * p.L_SALT / D_s)
    C = k_f_plus / (1.0 + k_f_minus * p.L_SALT / D_s)

    a = 2.0 * A
    b = C + D_m / p.L_NI
    c = -D_m * c_0 / p.L_NI

    c_metal = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)

    w_rec = A * c_metal**2
    w_f = C * c_metal
    flux = 2 * w_rec + w_f

    return c_metal, flux, w_f / flux


def build_model(
    k_r_plus,
    k_r_minus,
    k_f_plus,
    k_f_minus,
    temperature=p.TEMPERATURE,
    p_up=p.P_UP,
    final_time=FINAL_TIME,
):
    """Two-slab problem with the R and F channels sharing one interface."""
    metal = F.Material(D_0=p.D_0_NI, E_D=p.E_D_NI)
    salt = F.Material(
        D_0={"H2": p.D_0_SALT, "HF": p.D_0_SALT},
        E_D={"H2": p.E_D_SALT, "HF": p.E_D_SALT},
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

    H = F.Species("H", subdomains=[vol_metal])
    H2 = F.Species("H2", subdomains=[vol_salt])
    HF = F.Species("HF", subdomains=[vol_salt])
    model.species = [H, H2, HF]
    model.temperature = temperature

    # same id, so the two contributions add on the same facets: that is all
    # Model 3 requires
    model.interfaces = [
        F.InterfaceReaction(
            id=5,
            subdomains=[vol_metal, vol_salt],
            k_plus=k_r_plus,
            k_minus=k_r_minus,
            reactants=[H] * 2,
            products=[H2],
        ),
        F.InterfaceReaction(
            id=5,
            subdomains=[vol_metal, vol_salt],
            k_plus=k_f_plus,
            k_minus=k_f_minus,
            reactants=[H],
            products=[HF],
        ),
    ]

    model.boundary_conditions = [
        F.SievertsBC(
            subdomain=upstream, species=H, pressure=p_up,
            S_0=p.K_S_0_NI, E_S=p.E_K_S_NI,
        ),
        F.FixedConcentrationBC(subdomain=downstream, species=H2, value=0.0),
        F.FixedConcentrationBC(subdomain=downstream, species=HF, value=0.0),
    ]

    flux_h2 = F.SurfaceFlux(field=H2, surface=downstream)
    flux_hf = F.SurfaceFlux(field=HF, surface=downstream)
    profile_metal = F.Profile1DExport(field=H, subdomain=vol_metal)
    profile_h2 = F.Profile1DExport(field=H2, subdomain=vol_salt)
    profile_hf = F.Profile1DExport(field=HF, subdomain=vol_salt)
    model.exports = [flux_h2, flux_hf, profile_metal, profile_h2, profile_hf]

    model.settings = F.Settings(
        atol=1e10,
        rtol=1e-12,
        max_iterations=100,
        transient=True,
        final_time=final_time,
        stepsize=F.Stepsize(
            initial_value=INITIAL_STEP,
            growth_factor=1.2,
            cutback_factor=0.8,
            target_nb_iterations=4,
            max_stepsize=MAX_STEP,
        ),
    )

    return model, (flux_h2, flux_hf), (profile_metal, profile_h2, profile_hf)


def run(damkohler=DAMKOHLER, branching=1.0, **kwargs):
    """Solve one (Da, B) case.

    Returns a dict with the time histories and the steady-state diagnostics.
    The atomic flux is 2 J_H2 + J_HF, since each molecular carrier takes two
    hydrogen atoms out of the metal and each fluoride carrier takes one.
    """
    constants = channel_constants(damkohler, branching)
    model, fluxes, profiles = build_model(*constants, **kwargs)
    model.initialise()
    model.run()

    flux_h2, flux_hf = (np.abs(np.array(f.data)) for f in fluxes)
    profile_metal, profile_h2, profile_hf = profiles

    times = np.array(profile_metal.t)
    flux_atomic = 2.0 * flux_h2 + flux_hf

    c_metal = np.array([row[-1] for row in profile_metal.data])
    c_h2 = np.array([row[0] for row in profile_h2.data])
    c_hf = np.array([row[0] for row in profile_hf.data])

    # channel rates from the computed traces, i.e. from the same mass-action
    # expressions the solver assembled
    k_r_plus, k_r_minus, k_f_plus, k_f_minus = constants
    w_rec = k_r_plus * c_metal**2 - k_r_minus * c_h2
    w_f = k_f_plus * c_metal - k_f_minus * c_hf

    return {
        "times": times,
        "flux_atomic": flux_atomic,
        "flux_h2": flux_h2,
        "flux_hf": flux_hf,
        "hf_fraction": flux_hf / flux_atomic,
        "c_metal": c_metal,
        "w_rec": w_rec,
        "w_f": w_f,
        "branching_measured": w_f / (2.0 * w_rec),
        "constants": constants,
    }


def sweep(all_b, damkohler=DAMKOHLER):
    times_lte, flux_lte, _, _, _ = lte_baseline.run()
    lag_lte, steady_lte = time_lag(times_lte, flux_lte)

    print(f"Da = {damkohler:g} held fixed; LTE reference J_ss = {steady_lte:.4e}, "
          f"t_lag = {lag_lte:.4e} s\n")
    print(f"{'B target':>10s} {'B measured':>12s} {'J_ss/J_LTE':>11s} "
          f"{'HF fraction':>12s} {'t_lag/t_LTE':>12s}")

    results, rows = [], []
    for branching in all_b:
        out = run(damkohler, branching)
        lag, steady = time_lag(out["times"], out["flux_atomic"])

        results.append(out)
        rows.append([
            branching, out["branching_measured"][-1], steady, steady / steady_lte,
            out["hf_fraction"][-1], lag, lag / lag_lte, out["c_metal"][-1],
        ])

        print(f"{branching:10.2e} {out['branching_measured'][-1]:12.4e} "
              f"{steady / steady_lte:11.5f} {out['hf_fraction'][-1]:12.5f} "
              f"{lag / lag_lte:12.5f}")

    return (times_lte, flux_lte, steady_lte, lag_lte), results, np.array(rows)


def plot(reference, results, rows, filename):
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    times_lte, flux_lte, steady_lte, _ = reference
    b_measured = rows[:, 1]

    # Diverging about B = 1, where the two channels carry equal atomic flux, and
    # centred on grey rather than white so the middle curves do not disappear
    # into the background. Blue is recombination dominated, red is fluorination
    # dominated. A different family from the Da sweep keeps the two swept
    # variables visually distinct.
    cmap = LinearSegmentedColormap.from_list(
        "branching", ["#2166ac", "#8c8c8c", "#b2182b"]
    )
    # The measured branching ratio is not symmetric about 1: the F channel
    # drains the interfacial loading, so B climbs faster than the forward
    # constant does. A plain symmetric log scale would therefore leave most of
    # the blue half unused. Two slopes on log10(B) put grey exactly on B = 1
    # and still reach both ends of the map.
    log_b = np.log10(b_measured)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=log_b.min(), vmax=log_b.max())

    def colour(b):
        return cmap(norm(np.log10(b)))

    fig, (ax_shape, ax_split, ax_lag) = plt.subplots(
        1, 3, figsize=(11.5, 3.7), constrained_layout=True
    )

    window = 14.0  # hours

    # (a) each transient on its own steady value: the redox state alone retimes
    # the breakthrough, with every bulk property held fixed
    for b, out in zip(b_measured, results):
        steady = out["flux_atomic"][-1]
        ax_shape.plot(
            out["times"] / 3600, out["flux_atomic"] / steady,
            color=colour(b), linewidth=1.3,
        )
    ax_shape.plot(
        times_lte / 3600, flux_lte / steady_lte,
        color=GUIDE, linestyle="--", linewidth=1.4,
    )
    ax_shape.annotate(
        "LTE", xy=(11.0, 0.86), color=GUIDE, weight="bold", fontsize=9,
    )
    ax_shape.annotate(
        f"Da $= {DAMKOHLER:g}$ throughout",
        xy=(0.97, 0.06), xycoords="axes fraction", ha="right",
        fontsize=8, color="0.3",
    )
    ax_shape.set_xlabel("time [h]")
    ax_shape.set_ylabel(r"$J / J_\mathrm{steady}$")
    ax_shape.set_xlim(0, window)
    ax_shape.set_ylim(0, 1.05)
    ax_shape.set_title("(a) transient shape", loc="left", fontsize=10)

    # (b) the carrier split, a direct output of the model and measurable by
    # soluble/insoluble discrimination, which LTE cannot even define
    for b, out in zip(b_measured, results):
        ax_split.plot(
            out["times"] / 3600, out["hf_fraction"],
            color=colour(b), linewidth=1.3,
        )
    # every case starts HF-dominated and settles: B goes inversely with the
    # interfacial loading, Eq. (branching), and the loading builds up during the
    # transient. The split is a function of time, not a property of the pair
    ax_split.annotate(
        "split drifts as the\ninterfacial loading builds",
        xy=(1.1, 0.42), xytext=(4.6, 0.60),
        fontsize=8, color="0.3",
        arrowprops=dict(arrowstyle="->", color="0.5", linewidth=0.9,
                        connectionstyle="arc3,rad=0.25"),
    )
    ax_split.set_xlabel("time [h]")
    ax_split.set_ylabel("fraction of atomic flux carried by HF")
    ax_split.set_xlim(0, window)
    ax_split.set_ylim(0, 1.02)
    ax_split.set_title("(b) flux speciation", loc="left", fontsize=10)

    # (c) steady flux and time lag against the measured branching ratio
    ax_lag.semilogx(
        b_measured, rows[:, 3], color="C0", marker="o", markersize=4,
        markerfacecolor="none", linewidth=1.4, label=r"$J_\mathrm{ss}/J_\mathrm{LTE}$",
    )
    ax_lag.semilogx(
        b_measured, rows[:, 6], color="C2", marker="s", markersize=4,
        markerfacecolor="none", linewidth=1.4,
        label=r"$t_\mathrm{lag}/t_\mathrm{lag}^\mathrm{LTE}$",
    )
    ax_lag.axhline(1.0, color=GUIDE, linestyle="--", linewidth=1)
    ax_lag.set_xlabel(r"branching ratio $\mathcal{B}$ (measured)")
    ax_lag.set_title("(c) steady flux and lag", loc="left", fontsize=10)
    ax_lag.legend(frameon=False, loc="best", fontsize=8)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    bar = fig.colorbar(
        sm, ax=[ax_shape, ax_split],
        label=r"branching ratio $\mathcal{B}$ (measured)",
        location="bottom", pad=0.02, fraction=0.05, shrink=0.5, aspect=32,
    )
    # the bar carries log10(B), so label the decades as powers of ten
    decades = np.arange(np.ceil(log_b.min()), np.floor(log_b.max()) + 1)
    bar.set_ticks(decades)
    bar.set_ticklabels([f"$10^{{{int(d)}}}$" for d in decades])

    fig.savefig(filename)
    print(f"\nwrote {filename}")


if __name__ == "__main__":
    reference, results, rows = sweep(ALL_B)

    np.savetxt(
        "sweep_redox.csv",
        rows,
        delimiter=",",
        header="B_target,B_measured,J_steady,J_over_J_lte,hf_fraction,"
        "t_lag,t_lag_over_lte,c_metal_gamma",
        comments="",
    )
    plot(reference, results, rows, "sweep_redox.pdf")
