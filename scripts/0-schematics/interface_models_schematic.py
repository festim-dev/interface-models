"""The four interface models of Sec. "Interface conditions for hydrogen
transport", drawn on the same two-subdomain geometry so that only the interface
chemistry changes between panels.

  (a) LTE: nothing crosses. The two interfacial concentrations are tied by an
      algebraic constraint, so the interface has no rate constant, no timescale
      and no resistance, and the exponent of the constraint is chosen before the
      simulation is run.
  (b) Model 1, Eq. (model1_flux): one atom crosses without changing chemical
      identity. First order both ways.
  (c) Model 2, Eq. (model2_rate): two metal-side atoms recombine into one
      salt-side molecule. The 2:1 stoichiometry of Eq. (atom_balance) is drawn,
      the flux of atoms leaving the metal being twice the channel rate.
  (d) Model 3, Eqs. (model2_rate) + (model3_rate): the same atom has two exits
      of different order, and the branching ratio of Eq. (branching) sets how
      the atomic flux divides between them.

Purely illustrative: no solve is involved. Produces
interface_models_schematic.pdf, the roadmap figure of Sec. "A general kinetic
interface framework".
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

# the paper palette: metal side, molecular carrier, fluoride carrier
TEAL = "#1a4848"
AMBER = "#f7b000"
ORANGE = "#f46036"
GREEN = "#c9f2c7"

GAMMA = 5.0  # abscissa of the interface in every panel
WIDTH, HEIGHT = 10.0, 6.0
FLOOR = 1.6  # the diagram sits above this, the rate laws below it

# labels that land on an arrow or on the interface line need to punch through
LABEL_BOX = dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none")


def draw_domains(ax, left_label, right_label):
    """Two subdomains separated by the interface, with no axes of any kind."""
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, HEIGHT)
    ax.axis("off")

    ax.add_patch(
        Rectangle(
            (0, FLOOR), GAMMA, HEIGHT - FLOOR, facecolor=TEAL, alpha=0.11, linewidth=0
        )
    )
    ax.add_patch(
        Rectangle(
            (GAMMA, FLOOR),
            WIDTH - GAMMA,
            HEIGHT - FLOOR,
            facecolor=AMBER,
            alpha=0.17,
            linewidth=0,
        )
    )
    ax.plot([GAMMA, GAMMA], [FLOOR, HEIGHT], color=TEAL, linewidth=1.6, zorder=3)

    ax.text(0.35, HEIGHT - 0.45, left_label, color=TEAL, fontsize=9, va="center")
    ax.text(
        WIDTH - 0.35,
        HEIGHT - 0.45,
        right_label,
        color=TEAL,
        fontsize=9,
        va="center",
        ha="right",
    )
    ax.text(
        GAMMA,
        HEIGHT - 0.25,
        r"$\Gamma$",
        color=TEAL,
        fontsize=10,
        ha="center",
        va="center",
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none"),
    )


def title(ax, text):
    ax.text(0.15, HEIGHT + 0.35, text, color=TEAL, fontsize=9.5, weight="bold")


def caption(ax, lines):
    """Rate law and interface conditions, stacked under the diagram."""
    for i, line in enumerate(lines):
        ax.text(
            GAMMA,
            0.95 - 0.62 * i,
            line,
            color=TEAL,
            fontsize=8.5,
            ha="center",
            va="center",
        )


def atom(ax, x, y, label, color=TEAL, radius=0.4):
    ax.add_patch(
        Circle(
            (x, y),
            radius,
            facecolor=color,
            edgecolor="white",
            linewidth=1.0,
            zorder=5,
        )
    )
    ax.text(
        x,
        y,
        label,
        color="white",
        fontsize=8.5,
        ha="center",
        va="center",
        zorder=6,
        weight="bold",
    )


def diatomic(ax, x, y, left_label, right_label, colors=(AMBER, AMBER)):
    """Two bonded spheres, offset either side of x."""
    atom(ax, x - 0.3, y, left_label, color=colors[0], radius=0.38)
    atom(ax, x + 0.3, y, right_label, color=colors[1], radius=0.38)


def arrow(ax, start, end, color, rad=0.0, label=None, label_pos=None, style="-|>"):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(
            arrowstyle=style,
            color=color,
            linewidth=1.3,
            shrinkA=9,
            shrinkB=9,
            connectionstyle=f"arc3,rad={rad}",
        ),
        zorder=4,
    )
    if label is not None:
        ax.text(
            *label_pos,
            label,
            color=color,
            fontsize=8.5,
            ha="center",
            va="center",
            zorder=7,
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none"),
        )


def panel_lte(ax):
    """No arrows: that absence is the point of the panel."""
    draw_domains(ax, r"$\Omega_A$", r"$\Omega_B$")
    title(ax, "(a) local thermodynamic equilibrium")

    y = 3.4
    atom(ax, 3.3, y, "H")
    atom(ax, 6.7, y, "H", color=AMBER)

    # a rigid tie bar rather than a reaction: the constraint is instantaneous
    ax.plot([3.3, 6.7], [y, y], color=TEAL, linewidth=1.2, zorder=2)
    ax.add_patch(
        FancyBboxPatch(
            (GAMMA - 0.42, y - 0.32),
            0.84,
            0.64,
            boxstyle="round,pad=0.02",
            facecolor="white",
            edgecolor=TEAL,
            linewidth=1.2,
            zorder=6,
        )
    )
    ax.text(
        GAMMA, y, r"$=$", color=TEAL, fontsize=12, ha="center", va="center", zorder=7
    )

    caption(
        ax,
        [
            r"$c_A/K_{S,A} = c_B/K_{S,B}$",
            "no rate constant, no timescale, no resistance",
        ],
    )


def panel_model1(ax):
    draw_domains(ax, r"$\Omega_A$", r"$\Omega_B$")
    title(ax, "(b) Model 1: first-order exchange")

    y = 3.4
    atom(ax, 3.3, y, "H")
    atom(ax, 6.7, y, "H", color=AMBER)

    arrow(
        ax, (3.3, y), (6.7, y), TEAL, rad=-0.35, label=r"$k^{+}$", label_pos=(5.0, 4.35)
    )
    arrow(
        ax, (6.7, y), (3.3, y), AMBER, rad=-0.35, label=r"$k^{-}$", label_pos=(5.0, 2.20)
    )

    caption(
        ax,
        [
            r"$\phi = k^{+}c_A|_\Gamma - k^{-}c_B|_\Gamma$",
            r"$k^{+}/k^{-} = K_{S,B}/K_{S,A}$,  LTE as $k^{\pm}\to\infty$",
        ],
    )


def panel_model2(ax):
    draw_domains(ax, "metal", "salt")
    title(ax, "(c) Model 2: recombination into a carrier")

    atom(ax, 2.9, 4.7, "H")
    atom(ax, 2.9, 3.1, "H")
    diatomic(ax, 7.1, 4.0, "H", "H")

    arrow(ax, (2.9, 4.7), (6.7, 4.2), TEAL, rad=-0.10)
    arrow(ax, (2.9, 3.1), (6.7, 3.7), TEAL, rad=0.10)
    arrow(ax, (7.1, 3.4), (2.9, 3.1), AMBER, rad=0.30)

    ax.text(4.6, 4.95, r"$k_\mathrm{r}^{+}$", color=TEAL, fontsize=8.5, ha="center", zorder=7, bbox=LABEL_BOX)
    ax.text(5.0, 2.70, r"$k_\mathrm{r}^{-}$", color=AMBER, fontsize=8.5, ha="center", zorder=7, bbox=LABEL_BOX)

    caption(
        ax,
        [
            (
                r"$w_\mathrm{rec} = k_\mathrm{r}^{+}(c^\mathrm{m}_\mathrm{H})^2"
                r" - k_\mathrm{r}^{-}c^\mathrm{s}_{\mathrm{H_2}}$"
            ),
            (
                r"atoms leave at $2w_\mathrm{rec}$, one molecule enters at "
                r"$w_\mathrm{rec}$"
            ),
        ],
    )


def panel_model3(ax):
    draw_domains(ax, "metal", "salt")
    title(ax, "(d) Model 3: competing channels")

    # channel R, upper: second order, so it consumes a pair of metal-side atoms
    atom(ax, 2.8, 5.0, "H")
    atom(ax, 2.8, 3.4, "H")
    diatomic(ax, 7.2, 4.8, "H", "H")
    arrow(ax, (2.8, 5.0), (6.8, 4.9), TEAL, rad=-0.10)
    arrow(ax, (2.8, 3.4), (6.8, 4.5), TEAL, rad=-0.16)
    ax.text(4.4, 4.85, r"R", color=TEAL, fontsize=9, ha="center", weight="bold", zorder=7, bbox=LABEL_BOX)

    # channel F, lower: first order, and carries the salt redox state through a_F
    diatomic(ax, 7.2, 2.2, "H", "F", colors=(ORANGE, GREEN))
    arrow(ax, (2.8, 3.4), (6.8, 2.3), ORANGE, rad=0.14)
    ax.text(3.9, 2.50, r"F", color=ORANGE, fontsize=9, ha="center", weight="bold", zorder=7, bbox=LABEL_BOX)

    caption(
        ax,
        [
            (
                r"R: $2\,$H $\rightleftharpoons$ H$_2$   (quadratic)"
                r"     F: H $+\,a_\mathrm{F} \rightleftharpoons$ HF   (linear)"
            ),
            (
                r"$\mathcal{B} = w_\mathrm{F}/2w_\mathrm{rec}$ sets the split, "
                r"and it is kinetic"
            ),
        ],
    )


def plot(filename):
    mt.set_theme("urban")

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.4))
    for ax, panel in zip(
        axes.flat, (panel_lte, panel_model1, panel_model2, panel_model3)
    ):
        panel(ax)

    fig.subplots_adjust(
        left=0.01, right=0.99, top=0.94, bottom=0.02, wspace=0.06, hspace=0.30
    )
    fig.savefig(filename)


if __name__ == "__main__":
    plot("interface_models_schematic.pdf")
