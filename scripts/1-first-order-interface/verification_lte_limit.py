"""Verification of Model 1 against the closed-form two-slab solution, and of the
approach to LTE as the Damkoehler number grows.

Backs Sec. "Recovery of LTE in the fast-kinetics limit" and Sec. "Analytical
steady-state solutions" of the paper. Produces:

  - verification_lte_limit.csv : the raw table
  - verification_lte_limit.pdf : error against Da

Two quantities are reported per Da:

  err_analytical : relative difference between the FESTIM interfacial
                   concentrations and the closed form of analytical_solution.py.
                   This is a code-verification number and should sit at solver
                   tolerance (P1 elements are nodally exact for the piecewise
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

all_da = np.logspace(-2, 6, 17)

rows = []
for da in all_da:
    k_plus = da * D_A / L_A
    k_minus = alpha * k_plus

    model = run_model(D_A, D_B, k_plus, k_minus, x_interface, c_0, c_L, transient=False)

    c_A_num = model.exports[0].data[-1][-1]
    c_B_num = model.exports[1].data[-1][0]

    c_A_ana, c_B_ana, phi_ana = solution(
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
    ratio = c_A_num / c_B_num
    err_lte = abs(ratio - alpha) / alpha

    rows.append(
        (da, k_plus, c_A_num, c_B_num, c_A_ana, c_B_ana, err_analytical, err_lte)
    )
    print(
        f"Da={da:10.3e}  c_A={c_A_num:.8f} ({c_A_ana:.8f})  "
        f"c_B={c_B_num:.8f} ({c_B_ana:.8f})  "
        f"err_ana={err_analytical:.3e}  err_lte={err_lte:.3e}"
    )

rows = np.array(rows)

# observed order of convergence of err_lte in 1/Da between successive points
order = np.log(rows[:-1, 7] / rows[1:, 7]) / np.log(rows[1:, 0] / rows[:-1, 0])
print("\nobserved order in 1/Da (successive pairs):")
for da_lo, da_hi, p in zip(rows[:-1, 0], rows[1:, 0], order):
    print(f"  {da_lo:9.2e} -> {da_hi:9.2e} : p = {p:.4f}")

header = "Da,k_plus,c_A_num,c_B_num,c_A_ana,c_B_ana,err_analytical,err_lte"
np.savetxt(
    "verification_lte_limit.csv", rows, delimiter=",", header=header, comments=""
)

mt.set_theme("urban")
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
)
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

plt.figure(figsize=(6, 3.5))
plt.loglog(
    rows[:, 0], rows[:, 7], "o-", color="C0", markersize=4, label="departure from LTE"
)
plt.loglog(
    rows[:, 0],
    rows[0, 7] * rows[0, 0] / rows[:, 0],
    "--",
    color="C2",
    label=r"$\propto 1/\mathrm{Da}$",
)
plt.loglog(
    rows[:, 0], rows[:, 6], "s-", color="C1", markersize=4, label="error vs closed form"
)
plt.xlabel("Damköhler number (Da)")
plt.ylabel("relative error")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("verification_lte_limit.pdf")
