"""Verification of Model 1 against the closed-form two-slab solution, and of the
approach to LTE as the Damkoehler number grows.

Backs Sec. "Recovery of LTE in the fast-kinetics limit" and Sec. "Analytical
steady-state solutions" of the paper. Produces:

  - verification_lte_limit.csv : the raw table
  - verification_lte_limit.pdf : error against Da

Two quantities are reported per Da:

  err_analytical : relative difference between the FESTIM interfacial
                   concentrations and the closed form of analytical_solution.py.
                   This is a code-verification number and should sit at
                   round-off (P1 elements are nodally exact for the piecewise
                   linear steady profile).
  err_lte        : relative departure of the interfacial ratio c_A/c_B from the
                   LTE value k_minus/k_plus. This is a modelling number and
                   should decay as 1/Da.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np

from analytical_solution import solution
from example_usage import run_model

# same configuration as parametric_study.py
D_A = 0.5
D_B = 1.0
alpha = 0.5  # k_minus / k_plus
x_interface = 0.5
c_0 = 2.0
c_L = 1.0

L_A = x_interface
L_B = 1.0 - x_interface

COLUMNS = "Da,k_plus,c_A_num,c_B_num,c_A_ana,c_B_ana,err_analytical,err_lte"


def sweep(all_da):
    """Solve the steady two-slab problem at each Da and compare with the closed
    form. Returns one row per Da, with the columns of COLUMNS."""
    rows = []
    for da in all_da:
        k_plus = da * D_A / L_A
        k_minus = alpha * k_plus

        model = run_model(
            D_A, D_B, k_plus, k_minus, x_interface, c_0, c_L, transient=False
        )

        c_A_num = model.exports[0].data[-1][-1]
        c_B_num = model.exports[1].data[-1][0]

        c_A_ana, c_B_ana, _ = solution(
            c_0=c_0,
            c_L=c_L,
            D_A=D_A,
            D_B=D_B,
            L_A=L_A,
            L_B=L_B,
            k_plus=k_plus,
            k_minus=k_minus,
        )

        err_analytical = max(
            abs(c_A_num - c_A_ana) / abs(c_A_ana),
            abs(c_B_num - c_B_ana) / abs(c_B_ana),
        )
        err_lte = abs(c_A_num / c_B_num - alpha) / alpha

        rows.append(
            (da, k_plus, c_A_num, c_B_num, c_A_ana, c_B_ana, err_analytical, err_lte)
        )
        print(
            f"Da={da:10.3e}  c_A={c_A_num:.8f} ({c_A_ana:.8f})  "
            f"c_B={c_B_num:.8f} ({c_B_ana:.8f})  "
            f"err_ana={err_analytical:.3e}  err_lte={err_lte:.3e}"
        )

    return np.array(rows)


def plot(rows, filename):
    """Error against Da, labelled in the figure rather than in a legend box."""
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    da, err_ana, err_lte = rows[:, 0], rows[:, 6], rows[:, 7]

    plt.figure(figsize=(6, 3.5))

    # 1/Da guide, anchored on the last point so that it lies on the asymptote
    asymptotic = da >= 1
    plt.loglog(
        da[asymptotic],
        err_lte[-1] * da[-1] / da[asymptotic],
        linestyle="--",
        color="C1",
    )
    plt.annotate(
        r"$\propto 1/\mathrm{Da}$",
        xy=(da[-4], err_lte[-1] * da[-1] / da[-4]),
        xytext=(-52, -14),
        textcoords="offset points",
        color="C1",
        weight="bold",
    )

    plt.loglog(da, err_lte, marker="o", markersize=4, alpha=0.5, color="C0")
    plt.annotate(
        "departure from LTE",
        xy=(da[3], err_lte[3]),
        xytext=(6, 10),
        textcoords="offset points",
        color="C0",
        weight="bold",
    )

    plt.loglog(da, err_ana, marker="s", markersize=4, alpha=0.5, color="C2")
    plt.annotate(
        "error vs closed form",
        xy=(da[3], err_ana[3]),
        xytext=(6, 10),
        textcoords="offset points",
        color="C2",
        weight="bold",
    )

    plt.xlabel("Damköhler number (Da)")
    plt.ylabel("relative error")
    plt.tight_layout()
    plt.savefig(filename)


if __name__ == "__main__":
    rows = sweep(np.logspace(-2, 6, 17))

    # observed order of convergence of err_lte in 1/Da between successive points
    order = np.log(rows[:-1, 7] / rows[1:, 7]) / np.log(rows[1:, 0] / rows[:-1, 0])
    print("\nobserved order in 1/Da (successive pairs):")
    for da_lo, da_hi, p in zip(rows[:-1, 0], rows[1:, 0], order):
        print(f"  {da_lo:9.2e} -> {da_hi:9.2e} : p = {p:.4f}")

    np.savetxt(
        "verification_lte_limit.csv",
        rows,
        delimiter=",",
        header=COLUMNS,
        comments="",
    )
    plot(rows, "verification_lte_limit.pdf")
