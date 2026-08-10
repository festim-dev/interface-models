"""Why a large Damkoehler number does not by itself license LTE.

Two checks are usually run together on a kinetic interface: Da >> 1, and a
small defect ratio phi / (k_plus c_A|_Gamma). They are not the same check, and
this script exhibits a case where the first passes by four decades while the
second sits at unity.

The mechanism is visible in one line. Writing the Model 1 flux as a ratio of
the two directions,

    phi / (k_plus c_A|_Gamma) = 1 - c_B|_Gamma / (K c_A|_Gamma),

so the defect ratio is a function of the interfacial traces alone: it measures
how far the interface reaction sits from detailed balance. Da = k_plus L / D
carries the bulk lengths and diffusivities instead. Nothing ties them together,
and the boundary conditions can hold the trace ratio fixed while Da is swept.

From the analytical solution of analytical_solution.py, with u = c_L / (K c_0),
Da_A = k_plus L_A / D_A and Da_B = k_minus L_B / D_B,

    phi / (k_plus c_A|_Gamma) = (1 - u) / (1 + Da_B + u Da_A),

in which Da_A appears only multiplied by u. Sweep the downstream side clean
(u = 0) and Da_A drops out altogether, leaving 1 / (1 + Da_B). So a downstream
side that offers no back-pressure holds the defect ratio near unity however
large Da_A, and with it the two-sided group Da* = Da_A + Da_B, is made.

The script sweeps Da* over eight decades in two configurations that differ only
in the downstream diffusivity, both with the downstream face held at zero:

  balanced : D_B = D_A,      so Da_B = Da_A
  swept    : D_B = 1e6 D_A,  so Da_B = Da_A * 1e-6

and reports what LTE gets wrong in each. The flux error is 1 / (1 + Da*) in
both, as Eq. (model1_convergence) says it must be, and in the balanced
configuration the error on c_A|_Gamma is the same number. In the swept
configuration it is not: it is larger by 1e6, the diffusivity contrast, and at
Da* = 1e4 the analytical solutions differ by a factor of 101 on a flux the two agree on
to four digits.

The reading for Sec. "When can LTE still be used?": Da* bounds the error on the
flux, the defect ratio bounds the error on the interfacial state, and the
second check is redundant only if the flux is the sole output.

Produces:

  - defect_ratio_vs_damkohler.pdf

The analytical solution is the reference throughout, and FESTIM is run at three points
per configuration to confirm it, against both the kinetic InterfaceReaction and
the LTE Interface.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np

from analytical_solution import solution

# two slabs of equal thickness on [0, 1], upstream loaded, downstream swept
D_A = 1.0
x_interface = 0.5
L_A = x_interface
L_B = 1.0 - x_interface
c_0 = 1.0
c_L = 0.0

# K = k_plus / k_minus = K_S,B / K_S,A by Eq. (detailed_balance_1). Held at one
# so that the two configurations differ in D_B alone
K = 1.0

CONFIGURATIONS = {
    "balanced": 1.0 * D_A,
    "swept": 1e6 * D_A,
}

# Da at which the numbers are tabulated and the FESTIM checks are run
DA_HEADLINE = 1e4


def lte_solution(D_B):
    """Steady two-slab solution with an LTE interface, c_A|_Gamma = c_B|_Gamma / K.
    The interface is a pair of resistances in series with nothing between them,
    so this is the k_plus -> infinity limit of analytical_solution.solution."""
    R_A = L_A / D_A
    R_B = L_B / (K * D_B)
    phi = (c_0 - c_L / K) / (R_A + R_B)
    c_A_int = c_0 - phi * R_A
    return c_A_int, K * c_A_int, phi


def sweep(D_B, all_da):
    """Analytical solve at each Da, with Da the two-sided group Da* of
    Eq. (defect_1_explicit). Returns the defect ratio and the two errors LTE
    makes, on the flux and on the metal-side interfacial concentration."""
    R_A = L_A / D_A
    R_B = L_B / (K * D_B)
    c_A_lte, _, phi_lte = lte_solution(D_B)

    rows = []
    for da in all_da:
        # Da* = k_plus (L_A/D_A + L_B/(K D_B)) fixes k_plus, Eq. (defect_1_explicit)
        k_plus = da / (R_A + R_B)
        k_minus = k_plus / K

        c_A, c_B, phi = solution(c_0, c_L, D_A, D_B, L_A, L_B, k_plus, k_minus)

        rows.append(
            (
                da,
                k_plus,
                phi / (k_plus * c_A),  # the defect ratio, second check
                abs(phi - phi_lte) / phi_lte,  # error LTE makes on the flux
                abs(c_A - c_A_lte) / c_A_lte,  # error LTE makes on c_A|_Gamma
            )
        )

    return np.array(rows)


def festim_check(D_B, all_da):
    """Run both interface classes at each Da and compare with the analytical solutions.

    The LTE penalty term has to beat the diffusive stiffness of the stiffer
    slab, which is max(D)/h here, so it is scaled with D_B: the constraint
    error falls as 1/penalty_term and the default 10.0 leaves it uncontrolled
    once D_B/D_A is large.
    """
    import festim as F
    from dolfinx.log import set_log_level, LogLevel

    from example_usage import run_model

    R_A = L_A / D_A
    R_B = L_B / (K * D_B)
    c_A_lte, _, _ = lte_solution(D_B)

    for da in all_da:
        k_plus = da / (R_A + R_B)
        set_log_level(LogLevel.WARNING)
        model = run_model(
            D_A, D_B, k_plus, k_plus / K, x_interface, c_0, c_L, transient=False
        )
        c_A_num = model.exports[0].data[-1][-1]
        c_A_ana, _, _ = solution(c_0, c_L, D_A, D_B, L_A, L_B, k_plus, k_plus / K)

        c_A_lte_num = run_lte_model(D_B, penalty_term=1e6 * D_B)

        print(
            f"  Da*={da:8.1e}  kinetic c_A={c_A_num:.6e} (analytical solution "
            f"{c_A_ana:.6e}, rel {abs(c_A_num - c_A_ana) / c_A_ana:.1e})   "
            f"LTE c_A={c_A_lte_num:.6e} (analytical solution {c_A_lte:.6e}, rel "
            f"{abs(c_A_lte_num - c_A_lte) / c_A_lte:.1e})"
        )


def run_lte_model(D_B, penalty_term):
    """The same two-slab problem closed with F.Interface, the LTE baseline.
    Returns the metal-side interfacial concentration."""
    import festim as F

    mat_A = F.Material(D_0=D_A, E_D=0, K_S_0=1.0, E_K_S=0)
    mat_B = F.Material(D_0=D_B, E_D=0, K_S_0=K, E_K_S=0)

    vol1 = F.VolumeSubdomain1D(id=1, material=mat_A, borders=[0, x_interface])
    vol2 = F.VolumeSubdomain1D(id=2, material=mat_B, borders=[x_interface, 1])
    left = F.SurfaceSubdomain(id=3, locator=lambda x: np.isclose(x[0], 0))
    right = F.SurfaceSubdomain(id=4, locator=lambda x: np.isclose(x[0], 1))

    model = F.HydrogenTransportProblemDiscontinuous()
    model.subdomains = [vol1, vol2, left, right]
    model.mesh = F.Mesh1D(vertices=np.linspace(0, 1, 101))

    A = F.Species("A", subdomains=[vol1, vol2])
    model.species = [A]

    model.interfaces = [
        F.Interface(id=1, subdomains=[vol1, vol2], penalty_term=penalty_term)
    ]
    model.boundary_conditions = [
        F.FixedConcentrationBC(species=A, subdomain=left, value=c_0),
        F.FixedConcentrationBC(species=A, subdomain=right, value=c_L),
    ]
    model.temperature = 300
    model.settings = F.Settings(atol=1e-14, rtol=1e-14, transient=False)
    model.exports = [
        F.Profile1DExport(field=A, subdomain=vol1),
        F.Profile1DExport(field=A, subdomain=vol2),
    ]

    model.initialise()
    model.run()
    return model.exports[0].data[-1][-1]


def plot(results, filename):
    """Two panels on a shared Da* axis: the two checks, then the two errors."""
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    da = results["balanced"][:, 0]

    fig, (ax_check, ax_err) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

    # (a) the two checks of step 3. The Da check is 1/(1+Da*) by construction
    ax_check.loglog(da, 1 / (1 + da), linestyle="--", color="C1")
    ax_check.annotate(
        r"the $\mathrm{Da}$ check, $1/(1+\mathrm{Da}^{\star})$",
        xy=(da[14], 1 / (1 + da[14])),
        xytext=(-108, -20),
        textcoords="offset points",
        color="C1",
        weight="bold",
    )
    ax_check.loglog(
        da, results["balanced"][:, 2], marker="o", markersize=4, alpha=0.5, color="C0"
    )
    ax_check.annotate(
        "balanced",
        xy=(da[10], results["balanced"][10, 2]),
        xytext=(8, 12),
        textcoords="offset points",
        color="C0",
        weight="bold",
    )
    ax_check.loglog(
        da, results["swept"][:, 2], marker="s", markersize=4, alpha=0.5, color="C2"
    )
    ax_check.annotate(
        "swept downstream",
        xy=(da[10], results["swept"][10, 2]),
        xytext=(-30, 12),
        textcoords="offset points",
        color="C2",
        weight="bold",
    )
    ax_check.set_ylabel(r"$\phi\,/\,(k^{+}c_A|_\Gamma)$")

    # (b) what LTE actually gets wrong. The flux error is 1/(1+Da*) in both
    # configurations, and in the balanced one the interfacial concentration
    # carries the same error, so three of the four curves coincide
    ax_err.loglog(
        da, results["swept"][:, 3], marker="s", markersize=4, alpha=0.5, color="C0"
    )
    ax_err.loglog(
        da, results["balanced"][:, 3], marker="o", markersize=4, alpha=0.5, color="C0"
    )
    ax_err.loglog(
        da, results["balanced"][:, 4], marker="o", markersize=4, alpha=0.5, color="C0"
    )
    ax_err.annotate(
        "flux (both), and $c_A|_\\Gamma$ (balanced)",
        xy=(da[11], results["balanced"][11, 3]),
        xytext=(-46, -24),
        textcoords="offset points",
        color="C0",
        weight="bold",
    )
    ax_err.loglog(
        da, results["swept"][:, 4], marker="s", markersize=4, alpha=0.5, color="C2"
    )
    ax_err.annotate(
        "$c_A|_\\Gamma$, swept downstream",
        xy=(da[9], results["swept"][9, 4]),
        xytext=(-16, 14),
        textcoords="offset points",
        color="C2",
        weight="bold",
    )
    ax_err.set_xlabel(r"Damköhler number $\mathrm{Da}^{\star}$")
    ax_err.set_ylabel("relative error of the LTE condition")

    fig.tight_layout()
    fig.savefig(filename)


if __name__ == "__main__":
    all_da = np.logspace(-2, 6, 17)
    results = {name: sweep(D_B, all_da) for name, D_B in CONFIGURATIONS.items()}

    idx = int(np.argmin(abs(all_da - DA_HEADLINE)))
    print(f"at Da* = {all_da[idx]:.1e}\n")
    header = f"{'':26s}{'balanced':>14s}{'swept':>14s}"
    print(header)
    for label, column in [
        ("phi / (k+ c_A|Gamma)", 2),
        ("LTE error on flux", 3),
        ("LTE error on c_A|Gamma", 4),
    ]:
        print(
            f"  {label:24s}"
            f"{results['balanced'][idx, column]:>14.3e}"
            f"{results['swept'][idx, column]:>14.3e}"
        )

    print("\nFESTIM check against the analytical solutions:")
    for name, D_B in CONFIGURATIONS.items():
        print(f" {name}:")
        festim_check(D_B, [1e2, 1e4, 1e6])

    plot(results, "defect_ratio_vs_damkohler.pdf")
