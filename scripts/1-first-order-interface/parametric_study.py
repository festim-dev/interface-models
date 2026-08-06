"""Transient Damkoehler sweep of Model 1 on the two-slab problem.

Backs the transient figure of Sec. "Recovery of LTE in the fast-kinetics
limit": the interfacial ratio c_A/c_B against time, for Da from 1e-2 to 1e2,
approaching the LTE value k_minus/k_plus of Eq. (lte_ss). At small Da the
interface has a relaxation time of its own and the ratio spends the transient
far from that value, which is the timescale an algebraic closure denies it.
Produces:

  - parametric_study_damkohler.pdf

The solves come from `example_usage.run_model`, so the geometry, mesh and
species are exactly those of the steady sweep in `verification_lte_limit.py`
and the two figures can be read against each other.

This script previously also wrote parametric_study_damkohler_steady.pdf. That
panel is now the upper panel of verification_lte_limit.pdf, over a wider Da
range and from steady solves rather than the end of a transient, so only the
transient half is kept here.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np
from matplotlib.colors import LogNorm

from example_usage import run_model

# geometry and transport as in the steady sweep, so only the sweep differs
D_A = 0.5
D_B = 1.0
x_interface = 0.5
c_0 = 2.0
c_L = 1.0

# detailed balance, Eq. (detailed_balance_1): alpha = k_minus/k_plus = K_S,A/K_S,B
alpha = 0.5


def sweep(all_da):
    """Solve the transient problem at each Da and return the interfacial ratio
    against time. Da = k_plus L_A / D_A is an input here, the channel being
    first order (Sec. "Dimensionless criteria for LTE validity")."""
    histories = []
    for da in all_da:
        k_plus = da * D_A / x_interface
        model = run_model(D_A, D_B, k_plus, k_plus * alpha, x_interface, c_0, c_L)

        times = np.array(model.exports[0].t)
        ratios = np.array(
            [
                model.exports[0].data[idx][-1] / model.exports[1].data[idx][0]
                for idx in range(len(model.exports[0].data))
            ]
        )
        histories.append((times, ratios))
        print(f"Da={da:10.3e}  c_A/c_B at final time = {ratios[-1]:.8f}")

    return histories


def plot(all_da, histories, filename):
    """One curve per Da, coloured by Da against the colourbar, with the LTE
    value as a dashed line. Labelled in the figure rather than in a legend
    box, as elsewhere in these scripts."""
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    cmap = plt.get_cmap("RdYlGn")
    norm = LogNorm(vmin=min(all_da), vmax=max(all_da), clip=False)

    plt.figure()
    for da, (times, ratios) in zip(all_da, histories):
        plt.plot(times, ratios, color=cmap(norm(da)))

    plt.axhline(y=alpha, color="C0", linestyle="--")
    plt.annotate(
        f"$k_-/k_+ = {alpha}$",
        xy=(0.5, alpha),
        xytext=(0.3, alpha - 0.08),
        color="C0",
        weight="bold",
    )
    plt.annotate(
        r"Da $\ll 1$", xy=(0.1, 1.2), color=cmap(norm(min(all_da))), weight="bold"
    )
    plt.annotate(
        r"Da $\gg 1$",
        xy=(0.45, alpha + 0.1),
        color=cmap(norm(max(all_da))),
        weight="bold",
    )

    plt.xlabel("time")
    plt.ylabel("$c_A/c_B$")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # required for ScalarMappable
    plt.colorbar(sm, label="Damköhler number (Da)", ax=plt.gca())

    plt.tight_layout()
    plt.savefig(filename)


if __name__ == "__main__":
    all_da = np.logspace(-2, 2, 10)
    plot(all_da, sweep(all_da), "parametric_study_damkohler.pdf")
