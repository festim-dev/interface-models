"""The (Da, B) regime map for the validity of a local-equilibrium interface.

Purely analytical: no solve is involved. The two axes are the two ways an LTE
condition can fail, and they are independent.

  Da : Damkoehler number of a channel, Eq. (damkohler). A single channel is
       locally equilibrated when Da >> 1. The relative error made on the flux by
       assuming equilibrium is 1/(1 + Da), which is exact for Model 1 by
       Eq. (model1_convergence).

  B  : branching ratio between the two channels, Eq. (branching). A single
       algebraic law exists only when one channel dominates. The fraction of the
       atomic flux carried by the minority channel, min(1, B)/(1 + B), is the
       part of the flux that any single-channel model fails to describe.

The field plotted is the larger of the two, which is an indicator of how badly
the best available LTE condition does, not a rigorous error bound. The apparent
exponent n = (2 + B)/(1 + B) of Eq. (n_of_B) is carried on the right-hand axis.

Produces regime_map.pdf, the figure of Sec. "Dimensionless criteria for LTE
validity".
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import FormatStrFormatter, NullLocator


def n_of_B(B):
    """Apparent interfacial exponent, Eq. (n_of_B)."""
    return (2.0 + B) / (1.0 + B)


def B_of_n(n):
    """Inverse of n_of_B, used for the right-hand axis."""
    return (2.0 - n) / (n - 1.0)


def within_channel_error(Da):
    """Relative error on the flux from assuming a channel is equilibrated."""
    return 1.0 / (1.0 + Da)


def between_channel_error(B):
    """Fraction of the atomic flux carried by the minority channel."""
    return np.minimum(1.0, B) / (1.0 + B)


def plot(filename):
    mt.set_theme("urban")
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    da = np.logspace(-2, 4, 400)
    branching = np.logspace(-3, 3, 400)
    DA, BR = np.meshgrid(da, branching)

    error = np.maximum(within_channel_error(DA), between_channel_error(BR))

    # good (pale green) to bad (orange), taken from the paper palette
    cmap = LinearSegmentedColormap.from_list(
        "lte_validity", ["#c9f2c7", "#aceca1", "#f7b000", "#f46036"]
    )

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    levels = np.logspace(-3, np.log10(0.5), 25)
    field = ax.contourf(
        DA,
        BR,
        error,
        levels=levels,
        cmap=cmap,
        norm=LogNorm(vmin=1e-3, vmax=0.5),
        extend="both",
    )
    ax.contour(
        DA, BR, error, levels=[0.01, 0.1], colors="#1a4848", linewidths=0.8, alpha=0.6
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Damköhler number (Da)")
    ax.set_ylabel(r"branching ratio $\mathcal{B}$")

    labels = [
        ("kinetically\nlimited", 3e-2, 1.0),
        ("no algebraic\nlaw exists", 3e2, 1.0),
        ("Sieverts / Henry\n$n \\to 2$", 3e2, 3e-3),
        ("linear\n$n \\to 1$", 3e2, 3e2),
    ]
    for text, x, y in labels:
        ax.annotate(
            text,
            xy=(x, y),
            color="#1a4848",
            weight="bold",
            ha="center",
            va="center",
        )

    ax.axvline(1.0, color="#1a4848", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.axhline(1.0, color="#1a4848", linestyle=":", linewidth=0.8, alpha=0.5)

    # the exponent is a function of B alone, so it rides on the right-hand axis.
    # the secondary axis inherits the parent's log locator, hence the resets
    exponent_axis = ax.secondary_yaxis("right", functions=(n_of_B, B_of_n))
    exponent_axis.yaxis.set_minor_locator(NullLocator())
    exponent_axis.set_yticks([1.05, 1.25, 1.5, 1.75, 1.95])
    exponent_axis.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    exponent_axis.set_ylabel("apparent exponent $n$")

    bar = fig.colorbar(field, ax=ax, pad=0.16, ticks=[1e-3, 1e-2, 1e-1, 0.5])
    bar.ax.set_yticklabels(["0.1%", "1%", "10%", "50%"])
    bar.ax.minorticks_off()
    bar.set_label("LTE failure indicator")

    fig.tight_layout()
    fig.savefig(filename)


if __name__ == "__main__":
    plot("regime_map.pdf")
