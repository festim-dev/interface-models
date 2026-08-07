"""Experiment C: the measured (Da, B) regime map at the Ni/FLiBe operating point.

scripts/2-higher-order-reactions/regime_map.py draws the same two axes
analytically, from the two closed-form indicators

    within a channel   : 1 / (1 + Da)
    between channels   : min(1, B) / (1 + B)

and says plainly that the larger of the two is an indicator and not a rigorous
error bound. This script computes the error itself. Every point of the grid is a
full transient solve of Model 3 at the dimensional operating point of
parameters.py, read against the LTE baseline, so the question it answers is
whether the heuristic is a safe guide at a real operating point and where it is
not.

The indicator is deliberately *not* drawn on the figure. It measures a
different thing and overlaying it invites the reader to compare quantities that
are not comparable: the between-channel term is the flux fraction in the
minority channel, so it asks whether *some* single algebraic law could describe
the interface, whereas this map asks how wrong *the conventional choice* is, the
Sieverts/Henry closure of lte_baseline.py that a macroscopic code implements.
At large B those diverge, because almost all the flux is then fluorination and a
single linear law would do well while Sieverts/Henry is missing the dominant
pathway. `indicator` below is kept for the comparison in the working notes, and
at (Da, B) = (100, 100) it reports 1 % against a true error of 167 %.

Two errors are mapped, because an LTE closure can be wrong about the steady
state and about the timing independently:

    steady    : J_ss / J_ss^LTE - 1
    transient : t_lag / t_lag^LTE - 1

Both signed, not absolute. The two failure modes push the steady flux in
opposite directions, a slow interface holding it below LTE and a second channel
carrying flux LTE has no pathway for, so taking a modulus would fold the map
about a line that is itself the most interesting feature and would break the
zero contour into islands. Panel (a) is drawn from the closed form in
`sweep_redox.model3_steady_state`, which the 81 solves confirm to 8e-6, so it
needs no grid at all; panel (b) needs the transients and carries its sample
points.

Both axes are nominal groups evaluated at the upstream Sieverts loading, the same
convention parameters.py uses for Da. That keeps the grid rectangular and makes
the axes quantities a reader can estimate for their own system without solving
anything. The branching ratio measured from the converged rates runs higher
than its nominal value, by more than a decade at the top of the range, because
the fluorination channel drains the interfacial loading; the measured value is
carried in the CSV.

The upper axis converts Da into the forward rate constant it corresponds to at
this operating point, which is the number a rate measurement would report.

Produces:

  - regime_map_hyperion.csv
  - regime_map_hyperion.pdf
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

import parameters as p
import lte_baseline
import sweep_redox
from sweep_damkohler import time_lag

ALL_DA = np.logspace(-2, 2, 9)
ALL_B = np.logspace(-2, 2, 9)

CSV = "regime_map_hyperion.csv"
GUIDE = "0.45"

COLUMNS = (
    "Da,B_nominal,B_measured,J_steady,J_over_J_lte,t_lag,t_lag_over_lte,"
    "hf_fraction,c_metal_gamma"
)


def indicator(damkohler, branching):
    """The analytical indicator of Sec. "Dimensionless criteria for LTE validity".

    The larger of the within-channel and between-channel terms, which is what
    the Sec. 2.7 map plots.
    """
    return np.maximum(
        1.0 / (1.0 + damkohler), np.minimum(1.0, branching) / (1.0 + branching)
    )


def compute(all_da=ALL_DA, all_b=ALL_B):
    """One transient solve per grid point. Writes the CSV as it goes."""
    times_lte, flux_lte, _, _, _ = lte_baseline.run()
    lag_lte, steady_lte = time_lag(times_lte, flux_lte)
    print(f"LTE reference: J_ss = {steady_lte:.4e} m^-2 s^-1, "
          f"t_lag = {lag_lte:.4e} s")
    print(f"grid: {len(all_da)} x {len(all_b)} = {len(all_da) * len(all_b)} solves\n")

    rows = []
    for damkohler in all_da:
        for branching in all_b:
            out = sweep_redox.run(damkohler, branching)
            lag, steady = time_lag(out["times"], out["flux_atomic"])

            rows.append([
                damkohler, branching, out["branching_measured"][-1], steady,
                steady / steady_lte, lag, lag / lag_lte,
                out["hf_fraction"][-1], out["c_metal"][-1],
            ])

            print(f"Da={damkohler:8.3e}  B={branching:8.3e}  "
                  f"J/J_LTE={steady / steady_lte:8.4f}  "
                  f"t_lag/t_LTE={lag / lag_lte:7.4f}  "
                  f"HF={out['hf_fraction'][-1]:.4f}", flush=True)

            np.savetxt(CSV, np.array(rows), delimiter=",", header=COLUMNS,
                       comments="")

    return np.array(rows)


def steady_error_field(n=241):
    """Signed steady-flux error of the LTE closure, from the closed form.

    The steady state of Model 3 is algebraic, so this panel needs no solves and
    can be drawn as finely as we like. `sweep_redox.model3_steady_state` agrees
    with the 81 transient solves of the grid to 8e-6, so the map is verified
    everywhere the grid samples it.

    The error is kept *signed*. Its two failure modes push in opposite
    directions, a slow interface holding the flux below LTE and a second channel
    carrying flux LTE has no pathway for, so the zero contour is a real feature
    and not an artefact of taking a modulus.
    """
    da = np.logspace(-2, 2, n)
    b = np.logspace(-2, 2, n)
    da_grid, b_grid = np.meshgrid(da, b, indexing="ij")

    _, flux, _ = sweep_redox.model3_steady_state(da_grid, b_grid)
    return da_grid, b_grid, flux / p.lte_steady_state()[3] - 1.0


def plot(rows, filename, all_da=ALL_DA, all_b=ALL_B):
    mt.set_theme("urban")
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    shape = (len(all_da), len(all_b))
    da_grid = rows[:, 0].reshape(shape)
    b_grid = rows[:, 1].reshape(shape)
    lag_error = (rows[:, 6] - 1.0).reshape(shape)

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4), constrained_layout=True)

    # ---- (a) steady flux, signed, closed form on a fine grid ----------------
    da_fine, b_fine, steady = steady_error_field()

    # teal where LTE reads high, orange where it reads low, grey through zero so
    # the neutral band does not disappear into the page
    diverging = LinearSegmentedColormap.from_list(
        "lte_error", ["#1a4848", "#8c8c8c", "#f46036"]
    )
    # asymmetric on purpose: the error runs from about -0.9 to +1.7, so forcing
    # a symmetric range would leave half the teal unused and wash out the whole
    # interface-limited side of the map
    low, high = steady.min(), steady.max()
    mesh_a = axes[0].contourf(
        da_fine, b_fine, steady,
        levels=np.unique(
            np.concatenate([np.linspace(low, 0, 21), np.linspace(0, high, 21)])
        ),
        cmap=diverging, norm=TwoSlopeNorm(vcenter=0.0, vmin=low, vmax=high),
    )

    # the whole point of the panel: LTE is exact along this line, for the wrong
    # reason, and panel (b) shows the transient is not
    zero = axes[0].contour(
        da_fine, b_fine, steady, levels=[0.0], colors="white", linewidths=2.2
    )
    axes[0].clabel(zero, fmt={0.0: "LTE exact"}, fontsize=8.5)

    # naming the two sides is worth more than another pair of contours: the
    # failure modes push in opposite directions and that is why they can cancel
    for x, y, text, colour in (
        (0.028, 0.03, "LTE reads high\n(slow interface)", "#d8e6e4"),
        (2.2, 22.0, "LTE reads low\n(missing channel)", "#4a1a08"),
    ):
        axes[0].annotate(
            text, xy=(x, y), fontsize=8.5, weight="bold", color=colour,
            ha="left", va="center",
        )

    fig.colorbar(
        mesh_a, ax=axes[0], location="bottom", pad=0.02, fraction=0.06,
        shrink=0.85, aspect=26, ticks=[-0.75, -0.5, -0.25, 0.0, 0.5, 1.0, 1.5],
        label=r"$J_\mathrm{ss}/J_\mathrm{ss}^\mathrm{LTE} - 1$",
    )

    # ---- (b) time lag, from the 81 transient solves -------------------------
    mesh_b = axes[1].contourf(
        da_grid, b_grid, lag_error,
        levels=np.linspace(0.0, 0.5, 26), cmap="YlOrRd", extend="max",
    )
    lag_bands = axes[1].contour(
        da_grid, b_grid, lag_error, levels=[0.1], colors="#1a4848", linewidths=1.8
    )
    axes[1].clabel(lag_bands, fmt={0.1: "10 %"}, fontsize=8)
    axes[1].plot(
        rows[:, 0], rows[:, 1], linestyle="none", marker=".", markersize=2.0,
        color="0.35", alpha=0.55,
    )

    fig.colorbar(
        mesh_b, ax=axes[1], location="bottom", pad=0.02, fraction=0.06,
        shrink=0.85, aspect=26, ticks=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
        label=r"$t_\mathrm{lag}/t_\mathrm{lag}^\mathrm{LTE} - 1$",
    )

    for ax, title in (
        (axes[0], "(a) steady flux, closed form"),
        (axes[1], "(b) time lag, 81 transient solves"),
    ):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-2, 1e2)
        ax.set_ylim(1e-2, 1e2)
        ax.set_xlabel("Damköhler number Da")
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(False)

    axes[0].set_ylabel(r"branching ratio $\mathcal{B}$ (nominal)")
    axes[1].set_yticklabels([])

    # Da on the top axis as the rate constant it stands for at this operating
    # point. The recombination channel is second order in the metal-side
    # loading, so its forward constant is a m^4/s and not a velocity; the
    # fluorination channel's is a velocity. That the two channels' rate
    # constants do not even share units is the same point Sec. "A solubility
    # constant has no law-independent units" makes about solubilities.
    def da_to_k(damkohler):
        return p.k_plus_from_damkohler(np.asarray(damkohler))

    def k_to_da(k_plus):
        return p.damkohler_recombination(np.asarray(k_plus))

    for ax in axes:
        secondary = ax.secondary_xaxis("top", functions=(da_to_k, k_to_da))
        secondary.set_xlabel(r"$k_\mathrm{r}^+$ [m$^4$ s$^{-1}$]", fontsize=9)

    fig.savefig(filename)
    print(f"\nwrote {filename}")


if __name__ == "__main__":
    rows = compute()
    plot(rows, "regime_map_hyperion.pdf")
