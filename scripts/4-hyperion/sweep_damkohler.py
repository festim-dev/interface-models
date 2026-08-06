"""Experiment A: the permeation transient at fixed transport, varying Da.

Everything about the bulk is held still. The diffusivities, the solubilities,
the geometry and the upstream pressure are the ones in parameters.py and none
of them moves across this sweep. Only the rate constant of the single
recombination channel does, and detailed balance carries k_- along with it so
the thermodynamics is unchanged too.

What moves is the transient. Panel (a) shows the flux histories against the LTE
reference; panel (b) shows them rescaled by their own steady values, which
isolates the shape from the amplitude and is the point of the figure: a
permeation curve is not a readout of the diffusivities alone. Panel (c) turns
that into the number an experimentalist would quote, by running the classical
time-lag inversion on each synthetic transient and reporting the diffusivity it
returns. The model diffusivity is fixed, so every departure from unity there is
the interface being mistaken for bulk transport.

This is the dimensional counterpart of the sweeps in
scripts/1-first-order-interface/, which establish the same limit on a
dimensionless two-slab problem. Here the reader can see how large Da has to be
at a real operating point.

Produces:

  - sweep_damkohler.csv
  - sweep_damkohler.pdf
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np
from matplotlib.colors import LogNorm

import parameters as p
import lte_baseline
import kinetic_transient

# symmetric about Da = 1, so the colour scale can be centred on the crossover
# rather than on an arbitrary end of the sweep
ALL_DA = np.logspace(-2, 2, 9)

GUIDE = "0.45"

try:
    _trapezoid = np.trapezoid
except AttributeError:  # numpy < 2
    _trapezoid = np.trapz


def time_lag(times, flux, steady_fraction=0.999):
    """Classical time lag of a permeation transient.

    The cumulative permeated amount tends to Q(t) -> J_ss (t - t_lag), so once
    the flux is within `steady_fraction` of its final value the lag is read off
    as t - Q(t)/J_ss. Averaging over the converged tail rather than using a
    single point keeps the estimate insensitive to where the tail is judged to
    start.

    The absolute lag carries the first-order error of the backward Euler
    integration and is still drifting at coarse steps. Every lag reported here
    is divided by the LTE one, which converges far faster because the two runs
    share the integrator and the error largely cancels.

    Returns (t_lag, J_steady).
    """
    flux_steady = flux[-1]
    cumulative = np.concatenate([[0.0], np.cumsum(np.diff(times) * (flux[1:] + flux[:-1]) / 2)])

    converged = flux >= steady_fraction * flux_steady
    if not converged.any():
        return np.nan, flux_steady

    tail = converged & (times > 0)
    return float(np.mean(times[tail] - cumulative[tail] / flux_steady)), flux_steady


def apparent_diffusivity(t_lag):
    """Diffusivity a single-slab time-lag analysis would report.

    The textbook inversion for one slab of thickness L is D = L^2 / (6 t_lag).
    Applied to this two-layer stack it is not the diffusivity of anything, which
    is the point: it is what the standard analysis returns when handed the
    measured curve, and parameters.py never changes the real diffusivities.
    """
    return (p.L_NI + p.L_SALT) ** 2 / (6.0 * t_lag)


def sweep(all_da):
    """One transient per Damkoehler number, plus the LTE reference."""
    times_lte, flux_lte, _, _, _ = lte_baseline.run()
    lag_lte, steady_lte = time_lag(times_lte, flux_lte)
    print(f"LTE reference: J_ss = {steady_lte:.6e} m^-2 s^-1, "
          f"t_lag = {lag_lte:.6e} s, D_app = {apparent_diffusivity(lag_lte):.6e} m^2/s")
    print(f"{'Da':>10s} {'J_ss':>14s} {'J_ss/J_LTE':>11s} {'t_lag [s]':>12s} "
          f"{'t_lag/t_LTE':>12s} {'D_app/D_LTE':>12s}")

    histories, rows = [], []
    for damkohler in all_da:
        k_plus = p.k_plus_from_damkohler(damkohler)
        times, flux, c_metal, c_h2, _ = kinetic_transient.run(k_plus)
        lag, steady = time_lag(times, flux)

        histories.append((times, flux, steady))
        rows.append(
            [damkohler, k_plus, steady, steady / steady_lte, lag, lag / lag_lte,
             apparent_diffusivity(lag) / apparent_diffusivity(lag_lte),
             c_metal[-1], c_h2[-1]]
        )

        print(f"{damkohler:10.2e} {steady:14.6e} {steady / steady_lte:11.6f} "
              f"{lag:12.4e} {lag / lag_lte:12.4f} "
              f"{apparent_diffusivity(lag) / apparent_diffusivity(lag_lte):12.4f}")

    return (times_lte, flux_lte, steady_lte, lag_lte), histories, np.array(rows)


def plot(reference, histories, rows, filename):
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    times_lte, flux_lte, steady_lte, _ = reference
    all_da = rows[:, 0]

    # centred on Da = 1: the midpoint of a diverging map should sit on the
    # crossover between interface control and transport control, not halfway
    # along whatever range happened to be swept
    cmap = plt.get_cmap("RdYlGn")
    span = 10 ** np.ceil(np.abs(np.log10(all_da)).max())
    norm = LogNorm(vmin=1.0 / span, vmax=span)

    fig, (ax_flux, ax_shape, ax_lag) = plt.subplots(
        1, 3, figsize=(11.5, 3.7), constrained_layout=True
    )

    # the shape panel is zoomed to the rise; the tails are flat and identical
    # by construction, and on the full axis the curves sit on top of each other
    shape_window = 14.0  # hours

    # (a) the transients as they would be measured, against the LTE reference
    for damkohler, (times, flux, _) in zip(all_da, histories):
        ax_flux.plot(times / 3600, flux, color=cmap(norm(damkohler)), linewidth=1.3)
    ax_flux.plot(times_lte / 3600, flux_lte, color=GUIDE, linestyle="--", linewidth=1.4)
    ax_flux.annotate(
        "LTE",
        xy=(times_lte[-1] / 3600, steady_lte),
        xytext=(-6, 5),
        textcoords="offset points",
        ha="right",
        color=GUIDE,
        weight="bold",
        fontsize=9,
    )
    ax_flux.set_xlabel("time [h]")
    ax_flux.set_ylabel(r"downstream flux [m$^{-2}$ s$^{-1}$]")
    ax_flux.set_xlim(0, times_lte[-1] / 3600)
    ax_flux.set_ylim(bottom=0)
    ax_flux.set_title("(a) as measured", loc="left", fontsize=10)

    # (b) each curve on its own steady value, which leaves only the shape
    for damkohler, (times, flux, steady) in zip(all_da, histories):
        ax_shape.plot(
            times / 3600, flux / steady, color=cmap(norm(damkohler)), linewidth=1.3
        )
    ax_shape.plot(
        times_lte / 3600, flux_lte / steady_lte,
        color=GUIDE, linestyle="--", linewidth=1.4,
    )
    ax_shape.annotate(
        r"Da $\ll 1$",
        xy=(5.6, 0.62),
        color=cmap(norm(all_da.min())),
        weight="bold",
        fontsize=9,
    )
    ax_shape.annotate(
        r"Da $\gg 1$",
        xy=(1.5, 0.86),
        color=cmap(norm(all_da.max())),
        weight="bold",
        fontsize=9,
    )
    ax_shape.set_xlabel("time [h]")
    ax_shape.set_ylabel(r"$J / J_\mathrm{steady}$")
    ax_shape.set_xlim(0, shape_window)
    ax_shape.set_ylim(0, 1.05)
    ax_shape.set_title("(b) shape alone", loc="left", fontsize=10)

    # (c) what the standard inversion returns, against a fixed model diffusivity
    ax_lag.semilogx(all_da, rows[:, 5], color="C0", marker="o", markersize=4,
                    markerfacecolor="none", linewidth=1.4)
    ax_lag.axhline(1.0, color=GUIDE, linestyle="--", linewidth=1)
    ax_lag.annotate(
        "LTE",
        xy=(all_da.max(), 1.0),
        xytext=(-4, 6),
        textcoords="offset points",
        ha="right",
        color=GUIDE,
        weight="bold",
        fontsize=9,
    )
    ax_lag.set_xlabel("Damköhler number Da")
    ax_lag.set_ylabel(r"$t_\mathrm{lag} / t_\mathrm{lag}^\mathrm{LTE}$")
    ax_lag.set_title("(c) inferred time lag", loc="left", fontsize=10)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(
        sm, ax=[ax_flux, ax_shape], label="Damköhler number Da",
        location="bottom", pad=0.02, fraction=0.05, shrink=0.5, aspect=32,
    )

    fig.savefig(filename)
    print(f"\nwrote {filename}")


if __name__ == "__main__":
    reference, histories, rows = sweep(ALL_DA)

    np.savetxt(
        "sweep_damkohler.csv",
        rows,
        delimiter=",",
        header="Da,k_plus,J_steady,J_over_J_lte,t_lag,t_lag_over_lte,"
        "D_app_over_lte,c_metal_gamma,c_H2_gamma",
        comments="",
    )
    plot(reference, histories, rows, "sweep_damkohler.pdf")
