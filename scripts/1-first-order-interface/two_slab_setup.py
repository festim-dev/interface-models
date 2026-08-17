"""The two-slab verification problem of Sec. "Recovery of LTE in the
fast-kinetics limit", and the series-resistance reading of its analytical solution.

Left: the steady profile at three Damkoehler numbers, plotted on a common scale
(c_A on side A, c_B/K on side B) so that the LTE condition Eq. (lte_ss) is
continuity. The jump left at the interface is then exactly the defect

    Delta = c_A|G - c_B|G / K = phi / k+ = (c_0 - c_L/K) / (1 + Da*)

of Eq. (defect_1_explicit), which closes as Da grows but is nonzero at any
finite rate constant, because the interface is carrying a flux.

Right: the same analytical solution read as three resistances in series,
Eq. (model1_analytical). LTE is the short circuit of the middle one.

Analytical throughout, no solve: the profiles come from analytical_solution.py,
the same analytical solution that verification_lte_limit.py checks FESTIM against.
Produces two_slab_setup.pdf.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

from analytical_solution import solution

TEAL = "#1a4848"
AMBER = "#f7b000"
ORANGE = "#f46036"

# the configuration of verification_lte_limit.py and parametric_study.py
D_A, D_B = 0.5, 1.0
L_A, L_B = 0.5, 0.5
c_0, c_L = 2.0, 1.0
ALPHA = 0.5  # k_minus / k_plus, so K = k_plus / k_minus = 2
K = 1.0 / ALPHA

DA_SHOWN = (0.1, 1.0, 100.0)
DA_COLOURS = (ORANGE, AMBER, TEAL)


def profile(da):
    """Steady profile on the common scale, as two linear segments.

    Returns (x, u) with u = c_A on side A and u = c_B/K on side B, plus the
    interfacial jump.
    """
    k_plus = da * D_A / L_A
    c_A_int, c_B_int, phi = solution(
        c_0=c_0,
        c_L=c_L,
        D_A=D_A,
        D_B=D_B,
        L_A=L_A,
        L_B=L_B,
        k_plus=k_plus,
        k_minus=ALPHA * k_plus,
    )
    x = np.array([0.0, L_A, L_A, L_A + L_B])
    u = np.array([c_0, c_A_int, c_B_int / K, c_L / K])
    return x, u, c_A_int - c_B_int / K


def panel_profiles(ax):
    ax.add_patch(
        Rectangle((0, -1), L_A, 10, facecolor=TEAL, alpha=0.07, linewidth=0, zorder=0)
    )
    ax.add_patch(
        Rectangle(
            (L_A, -1), L_B, 10, facecolor=AMBER, alpha=0.11, linewidth=0, zorder=0
        )
    )
    ax.axvline(L_A, color=TEAL, linewidth=1.2, zorder=1)

    for da, colour in zip(DA_SHOWN, DA_COLOURS):
        x, u, _ = profile(da)
        # the two segments are drawn as one polyline with a break at the
        # interface, which is where the discontinuity belongs
        ax.plot(x[:2], u[:2], color=colour, linewidth=1.8, zorder=4)
        ax.plot(x[2:], u[2:], color=colour, linewidth=1.8, zorder=4)
        label_x = 0.33
        ax.text(
            label_x,
            c_0 - (c_0 - u[1]) * label_x / L_A + 0.07,
            f"$\\mathrm{{Da}} = {da:g}$",
            color=colour,
            fontsize=8,
            ha="center",
            va="bottom",
            zorder=6,
            bbox=dict(boxstyle="round,pad=0.1", facecolor="white", edgecolor="none"),
        )

    # the defect at the slowest channel, which is the quantity Eq. (defect_1)
    # calls Delta and the middle resistor of the network on the right
    _, u_slow, jump = profile(DA_SHOWN[0])
    ax.annotate(
        "",
        xy=(L_A + 0.035, u_slow[1]),
        xytext=(L_A + 0.035, u_slow[2]),
        arrowprops=dict(
            arrowstyle="<|-|>", color=ORANGE, linewidth=1.1, shrinkA=0, shrinkB=0
        ),
        zorder=5,
    )
    ax.text(
        L_A + 0.075,
        0.5 * (u_slow[1] + u_slow[2]),
        r"$\Delta = J/k^{+}$",
        color=ORANGE,
        fontsize=8.5,
        va="center",
        zorder=6,
    )

    ax.set_xlim(0, L_A + L_B)
    ax.set_ylim(0.2, 2.25)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$c_A$  and  $c_B/K$")
    ax.set_yticks([0.5, 1.0, 1.5, 2.0])
    ax.set_xticks([0, L_A, L_A + L_B])
    ax.set_xticklabels(["$0$", r"$\Gamma$", "$L$"])

    ax.text(0.03, 2.05, "$c_0$", color=TEAL, fontsize=8.5)
    ax.text(0.93, 0.56, "$c_L/K$", color=TEAL, fontsize=8.5)
    ax.text(
        0.05,
        0.33,
        r"$D_A = 0.5$, $D_B = 1$, $L_A = L_B = 0.5$, $K = 2$",
        color=TEAL,
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none"),
    )
    ax.set_title(
        "(a) steady profiles on a common scale", fontsize=9.5, loc="left", color=TEAL
    )


def resistor(ax, x, y, label, expression, colour, width=1.7, height=0.5):
    ax.add_patch(
        FancyBboxPatch(
            (x - 0.5 * width, y - 0.5 * height),
            width,
            height,
            boxstyle="round,pad=0.02",
            facecolor="white",
            edgecolor=colour,
            linewidth=1.4,
            zorder=4,
        )
    )
    ax.text(x, y, label, color=colour, fontsize=9, ha="center", va="center", zorder=5)
    ax.text(x, y - 0.55, expression, color=colour, fontsize=8, ha="center", va="center")


def panel_network(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    y = 3.6
    ax.plot([0.7, 9.3], [y, y], color=TEAL, linewidth=1.2, zorder=2)

    resistor(ax, 2.4, y, "$R_A$", r"$L_A/D_A$", TEAL)
    resistor(ax, 5.0, y, r"$R_\Gamma$", r"$1/k^{+}$", ORANGE)
    resistor(ax, 7.6, y, "$R_B/K$", r"$L_B/(K D_B)$", TEAL)

    for x, label in ((0.7, "$c_0$"), (9.3, "$c_L/K$")):
        ax.plot([x], [y], marker="o", markersize=5, color=TEAL, zorder=5)
        ax.text(x, y + 0.45, label, color=TEAL, fontsize=9, ha="center")

    # the flux is the current through the chain
    ax.annotate(
        "",
        xy=(1.55, y - 0.95),
        xytext=(0.75, y - 0.95),
        arrowprops=dict(arrowstyle="-|>", color=TEAL, linewidth=1.2),
    )
    ax.text(1.75, y - 0.95, r"$J$", color=TEAL, fontsize=9, va="center")

    # LTE short-circuits the interfacial resistance
    bypass = y + 1.0
    ax.plot(
        [4.0, 4.0, 6.0, 6.0],
        [y, bypass, bypass, y],
        color=ORANGE,
        linewidth=1.4,
        linestyle=(0, (4, 2)),
        zorder=6,
    )
    ax.text(
        5.0,
        bypass + 0.25,
        r"LTE: $R_\Gamma \to 0$",
        color=ORANGE,
        fontsize=8.5,
        ha="center",
    )

    ax.text(
        5.0,
        1.15,
        r"$J = \dfrac{c_0 - c_L/K}{R_A + R_\Gamma + R_B/K}$,"
        r"$\quad \mathrm{Da}^{\star} = (R_A + R_B/K)\,/\,R_\Gamma$",
        color=TEAL,
        fontsize=8.5,
        ha="center",
        va="center",
    )
    ax.set_title(
        "(b) the same solution as resistances in series",
        fontsize=9.5,
        loc="left",
        color=TEAL,
    )


def plot(filename):
    mt.set_theme("urban")
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    fig, (ax_left, ax_right) = plt.subplots(
        1, 2, figsize=(7.0, 3.1), width_ratios=[1.0, 1.15]
    )
    panel_profiles(ax_left)
    panel_network(ax_right)

    fig.subplots_adjust(left=0.08, right=0.99, top=0.88, bottom=0.14, wspace=0.18)
    fig.savefig(filename)


if __name__ == "__main__":
    plot("two_slab_setup.pdf")
