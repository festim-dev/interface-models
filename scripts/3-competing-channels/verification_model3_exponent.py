"""Verification of Model 3 on the two-slab problem: the apparent interfacial
exponent against the branching ratio.

Backs Eq. (n_of_B) and the right-hand axis of the (Da, B) regime map. Two
channels share one interface, recombination into H2 (quadratic in the metal-side
loading) and fluorination into HF (linear), so the total atomic flux leaving the
metal is neither a square law nor a linear one. Sweeping the upstream loading
moves the branching ratio B of Eq. (branching) across unity, and the exponent
measured from the flux should follow

    n = dln J / dln c_m|G = (2 + B) / (1 + B) .

Nothing about that is imposed on the solver: the two channels are declared as
two InterfaceReaction objects sharing an interface id, their residual
contributions add, and n is measured afterwards from the computed fluxes the way
a permeation experiment would measure it. Produces:

  - verification_model3_exponent.csv   : the raw table
  - verification_model3_exponent.pdf   : flux and exponent against the loading
  - verification_model3_branching.pdf  : both exponent readouts against B

Three quantities are measured at each loading:

  err_analytical : relative difference between the FESTIM interfacial
                   concentrations and the analytical solution of
                   analytical_solution_model3.py. A code-verification number;
                   the steady profiles are piecewise linear and so nodally exact
                   on P1, which leaves the two interface terms alone under test.
  n_flux         : logarithmic slope of the total atomic flux against the
                   interfacial loading, by centred differences along the sweep.
                   Compared with (2 + B)/(1 + B) evaluated from the computed
                   channel rates. The residual difference is finite-difference
                   truncation, second order in the sweep spacing.
  n_conc         : the same slope read off the salt-side interfacial inventory
                   2 c_H2|G + c_HF|G instead of the flux. Compared with
                   (2 + B r)/(1 + B r) with r = R_HF/R_H2, which is why the two
                   carriers are given deliberately different diffusivities here.

The salt side is swept, c_H2 = c_HF = 0 on the outer liquid face, which is the
regime the appendix derivation assumes. The reverse terms are active throughout
and are what makes this a check rather than a tautology: they renormalise the
forward constants without moving the exponent.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np

import festim as F

from analytical_solution_model3 import apparent_exponent, branching, solution

# geometry and metal transport as in the Model 1 and Model 2 tests, so the three
# sweeps are read on the same problem and only the interface changes
D_m = 0.5  # metal, atomic H
D_H2 = 1.0  # liquid, molecular carrier
D_HF = 0.25  # liquid, fluoride carrier; deliberately different, see n_conc
x_interface = 0.5

L_m = x_interface
L_s = 1.0 - x_interface

# swept salt: both carriers held at zero on the outer liquid face
c_L2 = 0.0
c_LF = 0.0

# recombination channel, with K = kr_plus/kr_minus = K_H/K_S**2 as in the
# Model 2 test (Eq. detailed_balance_2)
kr_plus = 1.0
kr_minus = 2.0

# fluorination channel. a_F is the sweep variable: it is an operating condition
# set by the salt redox state, not a material constant, and the implementation
# has no activity prefactor so it is folded into the forward constant by hand
kf_plus = 1.0
kf_minus = 2.0
ALL_A_F = (0.1, 1.0, 10.0)

# ratio of the two downstream resistances, which separates the flux readout from
# the concentration readout
R_RATIO = (L_s / D_HF) / (L_s / D_H2)

# renormalised forward constants of App. (analytical): with both carriers held at
# zero downstream, each channel rate is an explicit function of the metal-side
# trace alone, w_rec = A_REC c**2 and w_F = b_fluo(a_F) c. The recombination one
# carries no a_F, so the quadratic branch is common to the three sweeps
A_REC = kr_plus / (1 + kr_minus * L_s / D_H2)


def b_fluo(a_F):
    """Renormalised forward constant of the fluorination channel."""
    return kf_plus * a_F / (1 + kf_minus * L_s / D_HF)


COLUMNS = (
    "a_F,c_0,c_m,c_H2,c_HF,w_rec,w_F,J,B,"
    "n_flux,n_flux_pred,n_conc,n_conc_pred,err_analytical"
)


def run_model(c_0, a_F):
    """Steady two-slab problem with two channels on one interface.

    `reactants=[H] * 2` puts the square on the recombination rate and the
    stoichiometric factor two on the metal-side flux, Eq. (bc_metal). The two
    InterfaceReaction objects carry the same id, so their contributions add on
    the same facets; that is all Model 3 requires.
    """
    my_model = F.HydrogenTransportProblemDiscontinuous()

    metal = F.VolumeSubdomain1D(
        id=1, material=F.Material(D_0=D_m, E_D=0), borders=[0, x_interface]
    )
    # the two carriers move at different rates in the liquid, so the liquid
    # material carries a per-species diffusivity rather than a single one
    liquid = F.VolumeSubdomain1D(
        id=2,
        material=F.Material(D_0={"H2": D_H2, "HF": D_HF}, E_D={"H2": 0.0, "HF": 0.0}),
        borders=[x_interface, 1],
    )
    left = F.SurfaceSubdomain(id=3, locator=lambda x: np.isclose(x[0], 0))
    right = F.SurfaceSubdomain(id=4, locator=lambda x: np.isclose(x[0], 1))

    my_model.subdomains = [metal, liquid, left, right]
    my_model.mesh = F.Mesh1D(vertices=np.linspace(0, 1, 101))

    H = F.Species("H", subdomains=[metal])
    H2 = F.Species("H2", subdomains=[liquid])
    HF = F.Species("HF", subdomains=[liquid])
    my_model.species = [H, H2, HF]

    my_model.interfaces = [
        F.InterfaceReaction(  # R: 2 H(m) <-> H2(s)
            id=1,
            subdomains=[metal, liquid],
            k_plus=kr_plus,
            k_minus=kr_minus,
            reactants=[H] * 2,
            products=[H2],
        ),
        F.InterfaceReaction(  # F: H(m) <-> HF(s), a_F folded into k_plus
            id=1,
            subdomains=[metal, liquid],
            k_plus=kf_plus * a_F,
            k_minus=kf_minus,
            reactants=[H],
            products=[HF],
        ),
    ]

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(species=H, subdomain=left, value=c_0),
        F.FixedConcentrationBC(species=H2, subdomain=right, value=c_L2),
        F.FixedConcentrationBC(species=HF, subdomain=right, value=c_LF),
    ]

    my_model.temperature = 300
    my_model.settings = F.Settings(
        atol=1e-14, rtol=1e-14, transient=False, max_iterations=100
    )
    my_model.exports = [
        F.Profile1DExport(field=H, subdomain=metal),
        F.Profile1DExport(field=H2, subdomain=liquid),
        F.Profile1DExport(field=HF, subdomain=liquid),
    ]

    my_model.initialise()
    my_model.run()
    return my_model


def sweep(all_c_0, a_F):
    """Solve at each upstream loading and compare with the analytical solution. Returns
    one row per loading, with the columns of COLUMNS."""
    c_m, c_H2, c_HF, err = [], [], [], []

    for c_0 in all_c_0:
        model = run_model(c_0, a_F)

        num = (
            model.exports[0].data[-1][-1],
            model.exports[1].data[-1][0],
            model.exports[2].data[-1][0],
        )
        ana = solution(
            c_0=c_0,
            c_L2=c_L2,
            c_LF=c_LF,
            D_m=D_m,
            D_H2=D_H2,
            D_HF=D_HF,
            L_m=L_m,
            L_s=L_s,
            kr_plus=kr_plus,
            kr_minus=kr_minus,
            kf_plus=kf_plus,
            kf_minus=kf_minus,
            a_F=a_F,
        )[:3]

        c_m.append(num[0])
        c_H2.append(num[1])
        c_HF.append(num[2])
        err.append(max(abs(n / a - 1.0) for n, a in zip(num, ana)))

        print(
            f"a_F={a_F:5.2f}  c_0={c_0:9.3e}  c_m={num[0]:12.6e}  "
            f"c_H2={num[1]:12.6e}  c_HF={num[2]:12.6e}  err_ana={err[-1]:.3e}"
        )

    c_m = np.array(c_m)
    c_H2 = np.array(c_H2)
    c_HF = np.array(c_HF)

    # channel rates from the computed traces, i.e. from the same mass-action
    # expressions the solver assembled
    w_rec = kr_plus * c_m**2 - kr_minus * c_H2
    w_F = kf_plus * a_F * c_m - kf_minus * c_HF

    J = 2 * w_rec + w_F  # total atomic flux, Eq. (model3_total_flux)
    S = 2 * c_H2 + c_HF  # salt-side interfacial inventory
    B = branching(w_rec, w_F)

    # measured exponents: centred logarithmic slopes along the sweep
    log_c = np.log(c_m)
    n_flux = np.gradient(np.log(J), log_c)
    n_conc = np.gradient(np.log(S), log_c)

    return np.column_stack(
        [
            np.full_like(c_m, a_F),
            all_c_0,
            c_m,
            c_H2,
            c_HF,
            w_rec,
            w_F,
            J,
            B,
            n_flux,
            apparent_exponent(B),
            n_conc,
            apparent_exponent(B * R_RATIO),
            np.array(err),
        ]
    )


GUIDE = "0.45"  # neutral grey; C3/C4 of the palette are pale accents and
# disappear against the background when used for a line
TICKS = [1.0, 1.25, 1.5, 1.75, 2.0]


def set_style():
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False


def plot_loading(rows, filename):
    """The two panels read against the interfacial loading.

    (a) the total atomic flux against the interfacial loading, which is what the
    exponent below is the slope of: each sweep starts on its own linear branch
    and ends on the quadratic branch the three share. (b) the exponent against
    the same loading, one sweep per redox state, which is the shape a pressure
    sweep would report: lines are Eq. (n_of_B), markers the slopes measured from
    the FESTIM fluxes. The two share their x axis, so the bend in a flux curve
    sits above the rise of its exponent.
    """
    set_style()

    fig, (ax_flux, ax_load) = plt.subplots(
        2,
        1,
        sharex=True,
        height_ratios=[1.15, 1],
        figsize=(6, 7),
        layout="constrained",
    )

    # (a) the flux itself. The exponent of the panel below is the local slope of
    # these curves, so the two horizontal guides there are the two straight lines
    # here: the quadratic branch 2 w_rec, common to the three sweeps because a_F
    # does not enter the recombination channel, and a linear branch b_fluo(a_F) c
    # per sweep
    c_grid = np.logspace(np.log10(rows[:, 2].min()), np.log10(rows[:, 2].max()), 400)

    # the crossover of each sweep, where its two channels carry equal atomic flux
    c_cross = {a_F: b_fluo(a_F) / (2 * A_REC) for a_F in ALL_A_F}

    # drawn from a third of the earliest crossover upwards; below that it is far
    # under every sweep and only stretches the axis
    c_quad = c_grid[c_grid > min(c_cross.values()) / 3]
    ax_flux.loglog(
        c_quad,
        2 * A_REC * c_quad**2,
        color=GUIDE,
        linestyle="--",
        linewidth=1,
        label="recombination alone, $n = 2$",
    )

    for i, a_F in enumerate(ALL_A_F):
        sel = rows[:, 0] == a_F
        c_m, J = rows[sel, 2], rows[sel, 7]

        ax_flux.loglog(
            c_grid,
            2 * A_REC * c_grid**2 + b_fluo(a_F) * c_grid,
            color=f"C{i}",
            linewidth=1.2,
        )
        ax_flux.loglog(
            c_m[::4],
            J[::4],
            linestyle="none",
            marker="o",
            markersize=4,
            markerfacecolor="none",
            color=f"C{i}",
        )
        # the linear branch each sweep leaves, stopped a little past its own
        # crossover so the three guides do not overrun the panel
        c_lin = c_grid[c_grid < 3 * c_cross[a_F]]
        ax_flux.loglog(
            c_lin,
            b_fluo(a_F) * c_lin,
            color=GUIDE,
            linestyle=":",
            linewidth=1,
            label="fluorination alone, $n = 1$" if i == 0 else None,
        )

    ax_flux.set_ylabel(r"total atomic flux $J$")
    ax_flux.legend(
        loc="upper left", frameon=False, handletextpad=0.5, borderaxespad=0.2
    )

    # (b) exponent against loading, one sweep per redox state. The two exponents
    # LTE offers are the horizontal guides; every sweep spends most of its range
    # between them
    for level, label, x, ha in (
        (2.0, "Sieverts/Henry, $n = 2$", rows[:, 2].min(), "left"),
        (1.0, "linear, $n = 1$", rows[:, 2].max(), "right"),
    ):
        ax_load.axhline(y=level, color=GUIDE, linestyle="--", linewidth=1)
        ax_load.annotate(
            label,
            xy=(x, level),
            xytext=(8 if ha == "left" else -8, 6),
            textcoords="offset points",
            ha=ha,
            color=GUIDE,
            weight="bold",
        )

    for i, a_F in enumerate(ALL_A_F):
        sel = rows[:, 0] == a_F
        c_m, n_meas, n_pred = rows[sel, 2], rows[sel, 9], rows[sel, 10]

        ax_load.semilogx(c_m, n_pred, color=f"C{i}", linewidth=1.2)
        ax_load.semilogx(
            c_m[::4],
            n_meas[::4],
            linestyle="none",
            marker="o",
            markersize=4,
            markerfacecolor="none",
            color=f"C{i}",
        )

        # label each sweep where it crosses n = 3/2, which is where it is
        # steepest and furthest from its neighbours
        k = np.argmin(np.abs(n_pred - 1.5))
        ax_load.annotate(
            f"$a_F = {a_F:g}$",
            xy=(c_m[k], n_pred[k]),
            xytext=(8, 2),
            textcoords="offset points",
            color=f"C{i}",
            weight="bold",
        )

    ax_load.set_ylim(0.92, 2.08)
    ax_load.set_yticks(TICKS)
    ax_load.set_xlabel(r"interfacial loading $c^m_H|_\Gamma$")
    ax_load.set_ylabel("apparent exponent $n$")

    fig.savefig(filename)
    print(f"wrote {filename}")


def plot_branching(rows, filename):
    """The two readouts against the branching ratio itself.

    Both follow Eq. (n_of_B), the flux one in B and the inventory one in
    B R_HF/R_H2, so the curves are the same shape a constant factor apart on the
    log axis. Plotting the inventory readout against B R_HF/R_H2 instead would
    collapse the two onto one curve, which is prettier and reads as though the
    two measurements agree; they agree on the law and not on the number, so they
    are kept apart.
    """
    set_style()

    fig, ax_branching = plt.subplots(figsize=(6, 3.3), layout="constrained")

    b_grid = np.logspace(
        np.log10(rows[:, 8].min()) - 0.3, np.log10(rows[:, 8].max()) + 0.3, 400
    )
    ax_branching.semilogx(
        b_grid,
        apparent_exponent(b_grid),
        color="C0",
        linewidth=1.4,
        zorder=1,
        label=r"flux: $n = (2 + \mathcal{B})/(1 + \mathcal{B})$",
    )
    ax_branching.semilogx(
        b_grid,
        apparent_exponent(b_grid * R_RATIO),
        color="C2",
        linewidth=1.4,
        zorder=1,
        label=rf"inventory: $n = (2 + {R_RATIO:g}\mathcal{{B}})"
        rf"/(1 + {R_RATIO:g}\mathcal{{B}})$",
    )

    for column, marker in ((9, "o"), (11, "s")):
        ax_branching.semilogx(
            rows[::3, 8],
            rows[::3, column],
            linestyle="none",
            marker=marker,
            markersize=4,
            markerfacecolor="none",
            color="C0" if column == 9 else "C2",
            zorder=3,
        )

    # the two curves cross n = 3/2 a factor R_HF/R_H2 apart, which is the whole
    # difference between the two readouts
    ax_branching.annotate(
        "",
        xy=(1.0, 1.5),
        xytext=(1.0 / R_RATIO, 1.5),
        arrowprops=dict(arrowstyle="<->", color=GUIDE, linewidth=1.2),
    )
    ax_branching.annotate(
        rf"$R_{{HF}}/R_{{H_2}} = {R_RATIO:g}$",
        xy=(1.0, 1.5),
        xytext=(12, 5),
        textcoords="offset points",
        ha="left",
        color=GUIDE,
        weight="bold",
    )

    ax_branching.set_ylim(0.92, 2.08)
    ax_branching.set_yticks(TICKS)
    ax_branching.set_xlabel(r"branching ratio $\mathcal{B}$")
    ax_branching.set_ylabel("apparent exponent $n$")
    ax_branching.legend(
        loc="upper right", frameon=False, handletextpad=0.5, borderaxespad=0.2
    )

    fig.savefig(filename)
    print(f"wrote {filename}")


if __name__ == "__main__":
    all_c_0 = np.logspace(-3, 2, 61)
    rows = np.vstack([sweep(all_c_0, a_F) for a_F in ALL_A_F])

    # the centred differences are second order in the sweep spacing, so the
    # endpoints (one-sided in np.gradient) are excluded from the comparison
    interior = np.ones(len(rows), dtype=bool)
    for a_F in ALL_A_F:
        idx = np.flatnonzero(rows[:, 0] == a_F)
        interior[idx[[0, -1]]] = False

    d_flux = np.abs(rows[interior, 9] - rows[interior, 10])
    d_conc = np.abs(rows[interior, 11] - rows[interior, 12])

    print(f"\nsweep spacing: dln c_0 = {np.log(all_c_0[1] / all_c_0[0]):.4f}")
    print(f"max |err_analytical|            : {rows[:, 13].max():.3e}")
    print(f"max |n_flux - (2+B)/(1+B)|      : {d_flux.max():.3e}")
    print(f"max |n_conc - (2+Br)/(1+Br)|    : {d_conc.max():.3e}  (r = {R_RATIO:g})")
    print(
        f"n_flux range                    : {rows[:, 9].min():.4f} to {rows[:, 9].max():.4f}"
    )
    print(
        f"B range                         : {rows[:, 8].min():.3e} to {rows[:, 8].max():.3e}"
    )

    # n = 3/2 at B = 1 is the midpoint of Eq. (n_of_B); report the closest point
    k = np.argmin(np.abs(np.log(rows[:, 8])))
    print(
        f"closest point to B = 1          : B = {rows[k, 8]:.4f}, n = {rows[k, 9]:.4f}"
    )

    np.savetxt(
        "verification_model3_exponent.csv",
        rows,
        delimiter=",",
        header=COLUMNS,
        comments="",
    )
    plot_loading(rows, "verification_model3_exponent.pdf")
    plot_branching(rows, "verification_model3_branching.pdf")
