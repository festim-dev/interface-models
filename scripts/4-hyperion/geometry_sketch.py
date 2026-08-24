"""Geometry figure for Sec. "System description and model setup".

Two panels. (a) the HYPERION vessel as an axisymmetric half-section, drawn to
the coordinates of mesh.py in the HYPERION modelling repo. (b) the 1D reduction
the rest of Sec. 5 uses, which keeps the membrane and the pool and drops
everything else.

The point of putting them side by side is to make the reduction auditable. The
salt wets a sidewall as well as the membrane, and that sidewall carries a
permeation path in parallel with the one being modelled; the figure states the
two areas so the reader can size the omission. The multidimensional treatment
is in the companion analysis cited as `hyperion_multidim` in the paper. The 1D
model here is for demonstration, not for inferring properties from data, so the
omission is a caveat and not an error bar.

Produces:

  - geometry_sketch.pdf
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
from matplotlib.patches import Rectangle

from parameters import (
    L_NI,
    L_SALT,
    MEMBRANE_AREA,
    R_INNER,
    R_OUTER,
    SIDEWALL_AREA,
    VESSEL_LAYERS,
)

MM = 1e3  # metres to millimetres, the units the figure is labelled in

NICKEL = "#1a4848"
SALT = "#f7b000"
GAS = "#f5f5f2"
GUIDE = "0.45"


def draw_vessel(ax):
    """Axisymmetric half-section, radius to the right, axis at r = 0."""
    for name, y_bottom, y_top, kind in VESSEL_LAYERS:
        colour = {"nickel": NICKEL, "salt": SALT, "gas": GAS}[kind]
        ax.add_patch(
            Rectangle(
                (0, y_bottom * MM),
                R_INNER * MM,
                (y_top - y_bottom) * MM,
                facecolor=colour,
                edgecolor="none",
            )
        )

    # the outer wall runs the full height, and is the surface the 1D model drops
    y_total = VESSEL_LAYERS[-1][2]
    ax.add_patch(
        Rectangle(
            (R_INNER * MM, 0),
            (R_OUTER - R_INNER) * MM,
            y_total * MM,
            facecolor=NICKEL,
            edgecolor="none",
        )
    )

    ax.plot([0, 0], [0, y_total * MM], color=GUIDE, linestyle="-.", linewidth=1)

    # Labels sit clear of the vessel and are joined to their layer by a thin
    # leader, because the membrane and the pool are 2 and 5 mm thick on a 111 mm
    # vessel and their labels would otherwise land on top of each other.
    label_x = 50.0
    labels = (
        ("top cap", 0.1101, 110.0, NICKEL),
        ("cover gas", 0.069, 72.0, "0.25"),
        ("FLiBe pool", 0.02657, 46.0, "#8a6200"),
        ("membrane", 0.023, 30.0, NICKEL),
        ("upstream plenum", 0.012, 13.0, "0.25"),
        ("bottom cap", 0.001, 0.0, NICKEL),
    )
    for text, y_layer, y_label, colour in labels:
        ax.annotate(
            text,
            xy=(R_OUTER * MM, y_layer * MM),
            xytext=(label_x, y_label),
            textcoords="data",
            va="center",
            ha="left",
            fontsize=8,
            color=colour,
            arrowprops=dict(
                arrowstyle="-", color="0.55", linewidth=0.8, shrinkA=2, shrinkB=3
            ),
        )

    # bracket over the two layers the 1D model keeps, which is what panel (b) is.
    # It goes on the axis side so it does not compete with the layer labels.
    bracket_x = -7.0
    ax.plot(
        [bracket_x, bracket_x],
        [22.0, 29.14],
        color="#f46036",
        linewidth=2.5,
        solid_capstyle="butt",
    )
    for y in (22.0, 29.14):
        ax.plot([bracket_x, bracket_x + 2.5], [y, y], color="#f46036", linewidth=1.5)
    ax.text(
        bracket_x - 2.5,
        25.6,
        "kept in (b)",
        rotation=90,
        ha="right",
        va="center",
        fontsize=8,
        weight="bold",
        color="#f46036",
    )

    ax.set_xlim(-16, 74)
    ax.set_ylim(-6, 120)
    ax.set_aspect("equal")
    ax.set_xlabel("$r$ [mm]")
    ax.set_ylabel("$z$ [mm]")
    ax.set_xticks([0, 20, 39])
    ax.set_title("(a) HYPERION vessel", loc="left", fontsize=10)


def draw_reduction(ax):
    """The 1D two-slab model, drawn to the thicknesses the solver uses."""
    total = (L_NI + L_SALT) * MM

    ax.add_patch(Rectangle((0, 0), L_NI * MM, 1, facecolor=NICKEL, edgecolor="none"))
    ax.add_patch(
        Rectangle((L_NI * MM, 0), L_SALT * MM, 1, facecolor=SALT, edgecolor="none")
    )

    ax.text(
        L_NI * MM / 2,
        0.5,
        "Ni",
        ha="center",
        va="center",
        color="white",
        weight="bold",
    )
    ax.text(
        L_NI * MM + L_SALT * MM / 2,
        0.5,
        "FLiBe",
        ha="center",
        va="center",
        color="#4a3200",
        weight="bold",
    )

    ax.axvline(L_NI * MM, color="#f46036", linewidth=2.5, ymin=0.28, ymax=0.72)
    ax.annotate(
        r"$\Gamma$",
        xy=(L_NI * MM, -0.12),
        ha="center",
        va="top",
        fontsize=12,
        weight="bold",
        color="#f46036",
    )

    # the two gas-side boundary conditions, drawn as flux into and out of the stack
    for x, dx, text_x, ha, text in (
        (-0.75, 0.6, -2.2, "left", "H$_2$ at $p_\\mathrm{up}$"),
        (total + 0.15, 0.6, total + 2.2, "right", "swept, $p \\approx 0$"),
    ):
        ax.annotate(
            "",
            xy=(x + dx, 0.5),
            xytext=(x, 0.5),
            arrowprops=dict(arrowstyle="->", color="0.35", linewidth=1.4),
        )
        ax.text(text_x, 0.68, text, ha=ha, fontsize=8, color="0.25")

    # thicknesses above the stack, clear of everything else
    for x0, x1, text in (
        (0.0, L_NI * MM, f"{L_NI * MM:.3f} mm"),
        (L_NI * MM, total, f"{L_SALT * MM:.3f} mm"),
    ):
        ax.annotate(
            "",
            xy=(x0, 1.18),
            xytext=(x1, 1.18),
            arrowprops=dict(arrowstyle="<->", color=GUIDE, linewidth=1),
        )
        ax.text((x0 + x1) / 2, 1.26, text, ha="center", fontsize=8, color=GUIDE)

    ax.set_xlim(-2.3, total + 2.3)
    ax.set_ylim(-0.75, 1.55)
    ax.set_yticks([])
    ax.set_xlabel("$x$ [mm]")
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.set_title("(b) 1D reduction", loc="left", fontsize=10)

    # ax.text(
    #     total / 2,
    #     -0.62,
    #     f"sidewall dropped: {SIDEWALL_AREA:.2e} m$^2$ wetted, against "
    #     f"{MEMBRANE_AREA:.2e} m$^2$ of membrane",
    #     ha="center",
    #     fontsize=7.5,
    #     color=GUIDE,
    # )


def plot(filename):
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    fig, (ax_vessel, ax_1d) = plt.subplots(
        1, 2, figsize=(8.6, 4.4), gridspec_kw={"width_ratios": [1, 1.15]}
    )
    draw_vessel(ax_vessel)
    draw_reduction(ax_1d)

    fig.tight_layout()
    fig.savefig(filename)
    print(f"wrote {filename}")
    print(f"  sidewall / membrane area = {SIDEWALL_AREA / MEMBRANE_AREA:.3f}")
    print(
        f"  sidewall / wetted Ni     = "
        f"{SIDEWALL_AREA / (SIDEWALL_AREA + MEMBRANE_AREA):.3f}"
    )


if __name__ == "__main__":
    plot("geometry_sketch.pdf")
