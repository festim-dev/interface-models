"""Verification of Model 2 on a metal/liquid problem: the same eight-decade
Damkoehler sweep as the Model 1 test, converging to the Sieverts/Henry condition.

Backs Sec. "Recovery of LTE in the fast-kinetics limit". The geometry is the
two-slab problem of `1-first-order-interface/`, but the interface now carries
the recombination channel 2 H(m) <-> H2(s) of Eq. (model2_rate): the metal
holds an atomic species, the liquid a molecular one, and the two are different
species rather than the same species on two sides. Produces:

  - verification_model2_lte_limit.csv : the raw table
  - verification_model2_lte_limit.pdf : the figure

The channel is second order, so the exchange velocity
2 k_plus c_m|G depends on the solution and Da is not an input to the problem.
Following Sec. "Dimensionless criteria for LTE validity" the sweep is
controlled by the upstream loading,

    Da = 2 k_plus c_0 L_m / D_m

with c_0 the imposed metal-side concentration standing in for the Sieverts
value K_S sqrt(P_up), and the value actually attained at the interface,

    Da_local = 2 k_plus c_m|G L_m / D_m

is reported as a diagnostic. Three quantities are measured per Da:

  err_analytical : relative difference between the FESTIM interfacial
                   concentrations and the analytical solution of
                   analytical_solution_model2.py. A code-verification number;
                   the steady profiles are piecewise linear and so nodally
                   exact on P1, which leaves the interface term alone under
                   test.
  err_lte        : relative departure of c_s|G from K c_m|G**2, i.e. from the
                   Sieverts/Henry condition Eq. (lte_sh). A modelling number,
                   equal to w/(k_plus c_m|G**2) by Eq. (defect_2), so it should
                   decay as 1/Da.
  da_local       : the diagnostic above, as a fraction of the control Da.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np

import festim as F

from analytical_solution_model2 import lte_solution, solution

# geometry and transport as in the Model 1 test, so the two sweeps are read on
# the same problem and only the interface channel differs
D_m = 0.5  # metal, atomic H
D_s = 1.0  # liquid, molecular H2
x_interface = 0.5
c_0 = 2.0  # atomic concentration imposed on the outer metal face
c_L = 1.0  # molecular concentration imposed on the outer liquid face

# detailed balance, Eq. (detailed_balance_2): K = k_plus/k_minus = K_H/K_S**2
K = 0.5

L_m = x_interface
L_s = 1.0 - x_interface

# loading used to convert the control Da into a rate constant
DA_SCALE = 2 * c_0 * L_m / D_m

COLUMNS = (
    "Da,Da_local,k_plus,c_m_num,c_s_num,c_m_ana,c_s_ana,err_analytical,err_lte"
)


def run_model(k_plus, k_minus):
    """Steady two-slab problem with a single recombination channel at the
    interface. `reactants=[H] * 2` is what puts the stoichiometric factor two
    on the metal-side flux and the square on the rate."""
    my_model = F.HydrogenTransportProblemDiscontinuous()

    metal = F.VolumeSubdomain1D(
        id=1, material=F.Material(D_0=D_m, E_D=0), borders=[0, x_interface]
    )
    liquid = F.VolumeSubdomain1D(
        id=2, material=F.Material(D_0=D_s, E_D=0), borders=[x_interface, 1]
    )
    left = F.SurfaceSubdomain(id=3, locator=lambda x: np.isclose(x[0], 0))
    right = F.SurfaceSubdomain(id=4, locator=lambda x: np.isclose(x[0], 1))

    my_model.subdomains = [metal, liquid, left, right]
    my_model.mesh = F.Mesh1D(vertices=np.linspace(0, 1, 101))

    H = F.Species("H", subdomains=[metal])
    H2 = F.Species("H2", subdomains=[liquid])
    my_model.species = [H, H2]

    my_model.interfaces = [
        F.InterfaceReaction(
            id=1,
            subdomains=[metal, liquid],
            k_plus=k_plus,
            k_minus=k_minus,
            reactants=[H] * 2,
            products=[H2],
        ),
    ]

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(species=H, subdomain=left, value=c_0),
        F.FixedConcentrationBC(species=H2, subdomain=right, value=c_L),
    ]

    my_model.temperature = 300
    my_model.settings = F.Settings(
        atol=1e-12, rtol=1e-12, transient=False, max_iterations=100
    )
    my_model.exports = [
        F.Profile1DExport(field=H, subdomain=metal),
        F.Profile1DExport(field=H2, subdomain=liquid),
    ]

    my_model.initialise()
    my_model.run()
    return my_model


def sweep(all_da):
    """Solve at each control Da and compare with the analytical solution. One row per
    Da, with the columns of COLUMNS."""
    rows = []
    for da in all_da:
        k_plus = da / DA_SCALE
        k_minus = k_plus / K

        model = run_model(k_plus, k_minus)

        c_m_num = model.exports[0].data[-1][-1]
        c_s_num = model.exports[1].data[-1][0]

        c_m_ana, c_s_ana, _ = solution(
            c_0=c_0,
            c_L=c_L,
            D_m=D_m,
            D_s=D_s,
            L_m=L_m,
            L_s=L_s,
            k_plus=k_plus,
            k_minus=k_minus,
        )

        err_analytical = max(
            abs(c_m_num - c_m_ana) / abs(c_m_ana),
            abs(c_s_num - c_s_ana) / abs(c_s_ana),
        )
        err_lte = abs(c_s_num / (K * c_m_num**2) - 1.0)
        da_local = 2 * k_plus * c_m_num * L_m / D_m

        rows.append(
            (
                da,
                da_local,
                k_plus,
                c_m_num,
                c_s_num,
                c_m_ana,
                c_s_ana,
                err_analytical,
                err_lte,
            )
        )
        print(
            f"Da={da:10.3e}  Da_loc={da_local:10.3e}  "
            f"c_m={c_m_num:.8f} ({c_m_ana:.8f})  "
            f"c_s={c_s_num:.8f} ({c_s_ana:.8f})  "
            f"err_ana={err_analytical:.3e}  err_lte={err_lte:.3e}"
        )

    return np.array(rows)


def plot(rows, filename):
    """Two panels on a shared Da axis, labelled in the figure rather than in a
    legend box: (a) the Sieverts/Henry readout approaching its LTE value,
    (b) the rate at which it gets there. Deliberately the same layout as the
    Model 1 figure, so the two sweeps can be read side by side.

    da_local is not plotted. It is a modelling diagnostic rather than a
    verification result, and it stays in the CSV and in the printed output.
    """
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    da = rows[:, 0]
    c_m, c_s = rows[:, 3], rows[:, 4]
    err_ana, err_lte = rows[:, 7], rows[:, 8]
    readout = c_s / c_m**2

    fig, (ax_ratio, ax_err) = plt.subplots(
        2, 1, figsize=(6, 5), sharex=True, height_ratios=[1, 1.3]
    )

    # (a) the physical readout: c_s/c_m^2 -> K is Eq. (lte_sh)
    ax_ratio.semilogx(da, readout, marker="o", markersize=4, alpha=0.5, color="C0")
    ax_ratio.axhline(y=K, color="C1", linestyle="--")
    ax_ratio.annotate(
        "$K_H/K_S^2$",
        xy=(da[0], K),
        xytext=(0, -16),
        textcoords="offset points",
        color="C1",
        weight="bold",
    )
    ax_ratio.annotate(
        "kinetically limited",
        xy=(da[1], readout[1]),
        xytext=(6, -14),
        textcoords="offset points",
        color="C0",
        weight="bold",
    )
    ax_ratio.annotate(
        "Sieverts/Henry recovered",
        xy=(da[-5], readout[-5]),
        xytext=(-90, 12),
        textcoords="offset points",
        color="C0",
        weight="bold",
    )
    ax_ratio.set_ylim(0, 0.62)
    ax_ratio.set_ylabel("$c^s_{H_2}/(c^m_H)^2$")

    # (b) the convergence rate. 1/Da guide anchored on the last point so that
    # it lies on the asymptote rather than on the pre-asymptotic first point
    asymptotic = da >= 1
    ax_err.loglog(
        da[asymptotic],
        err_lte[-1] * da[-1] / da[asymptotic],
        linestyle="--",
        color="C1",
    )
    ax_err.annotate(
        r"$\propto 1/\mathrm{Da}$",
        xy=(da[-4], err_lte[-1] * da[-1] / da[-4]),
        xytext=(-52, -14),
        textcoords="offset points",
        color="C1",
        weight="bold",
    )

    ax_err.loglog(da, err_lte, marker="o", markersize=4, alpha=0.5, color="C0")
    ax_err.annotate(
        "departure from Sieverts/Henry",
        xy=(da[3], err_lte[3]),
        xytext=(6, 10),
        textcoords="offset points",
        color="C0",
        weight="bold",
    )

    # several points are at exactly zero, which a log axis drops; floor them at
    # the machine epsilon they are indistinguishable from anyway
    ax_err.loglog(
        da,
        np.maximum(err_ana, np.finfo(float).eps),
        marker="s",
        markersize=4,
        alpha=0.5,
        color="C2",
    )
    ax_err.annotate(
        "error vs analytical",
        xy=(da[8], np.finfo(float).eps),
        xytext=(-30, 10),
        textcoords="offset points",
        color="C2",
        weight="bold",
    )

    ax_err.set_xlabel("Damköhler number (Da)")
    ax_err.set_ylabel("relative error")

    fig.tight_layout()
    fig.savefig(filename)


if __name__ == "__main__":
    rows = sweep(np.logspace(-2, 6, 17))

    # observed order of convergence of err_lte in 1/Da between successive points
    order = np.log(rows[:-1, 8] / rows[1:, 8]) / np.log(rows[1:, 0] / rows[:-1, 0])
    print("\nobserved order in 1/Da (successive pairs):")
    for da_lo, da_hi, p in zip(rows[:-1, 0], rows[1:, 0], order):
        print(f"  {da_lo:9.2e} -> {da_hi:9.2e} : p = {p:.4f}")

    c_m_lte, c_s_lte, w_lte = lte_solution(
        c_0=c_0, c_L=c_L, D_m=D_m, D_s=D_s, L_m=L_m, L_s=L_s, K=K
    )
    print(
        f"\nSieverts/Henry limit: c_m={c_m_lte:.8f}  c_s={c_s_lte:.8f}  w={w_lte:.8f}"
    )
    print(f"  approached to {abs(rows[-1, 3] / c_m_lte - 1):.3e} (c_m) at Da={rows[-1, 0]:.0e}")
    # asymptote of err_lte from Eq. (defect_2): w/(k_plus c_m^2), with w and
    # c_m at their LTE values
    print(f"  predicted err_lte asymptote: {w_lte * DA_SCALE / c_m_lte**2:.6f}/Da")

    np.savetxt(
        "verification_model2_lte_limit.csv",
        rows,
        delimiter=",",
        header=COLUMNS,
        comments="",
    )
    plot(rows, "verification_model2_lte_limit.pdf")
