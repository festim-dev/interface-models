"""Zero-flux test of the isotopologue channels.

Two mobile atomic species (H, T) in a metal feed three molecular carriers
(H2, HT, T2) in a liquid through three InterfaceReaction channels. The metal
side is held at fixed concentrations and the liquid has no outlet, so at steady
state every channel rate vanishes and the interface sits at true equilibrium.

We then measure the isotope exchange constant

    K_exch = c_HT^2 / (c_H2 * c_T2)

which the classical (mass-independent) statistical result puts at 4, i.e.
H2 : HT : T2 = 1 : 2 : 1 for equal H and T populations.

The point of the test is to find out what FESTIM's mass-action convention
gives, rather than to assume it.
"""

import matplotlib

matplotlib.use("Agg")

import festim as F
import matplotlib.pyplot as plt
import morethemes as mt
import numpy as np


def run(k_plus_HT, k_plus_homo=1.0, k_minus=1.0, c_H=1.0, c_T=1.0, duplicate_HT=False):
    """Equilibrate the three channels and return (c_H2, c_HT, c_T2) at Gamma.

    If duplicate_HT is True, the mixed channel is declared twice, as H + T and
    as T + H, to test whether declaring both orderings reproduces the
    statistical factor without touching the rate constants.
    """
    model = F.HydrogenTransportProblemDiscontinuous()

    metal = F.Material(D_0=1.0, E_D=0)
    liquid = F.Material(D_0=1.0, E_D=0)

    vol_m = F.VolumeSubdomain1D(id=1, material=metal, borders=[0, 0.5])
    vol_l = F.VolumeSubdomain1D(id=2, material=liquid, borders=[0.5, 1])
    left = F.SurfaceSubdomain(id=3, locator=lambda x: np.isclose(x[0], 0))

    model.subdomains = [vol_m, vol_l, left]
    model.mesh = F.Mesh1D(vertices=np.linspace(0, 1, 101))

    H = F.Species("H", subdomains=[vol_m])
    T = F.Species("T", subdomains=[vol_m])
    H2 = F.Species("H2", subdomains=[vol_l])
    HT = F.Species("HT", subdomains=[vol_l])
    T2 = F.Species("T2", subdomains=[vol_l])
    model.species = [H, T, H2, HT, T2]

    # 2H <-> H2, H + T <-> HT, 2T <-> T2
    model.interfaces = [
        F.InterfaceReaction(
            id=1,
            subdomains=[vol_m, vol_l],
            k_plus=k_plus_homo,
            k_minus=k_minus,
            reactants=[H, H],
            products=[H2],
        ),
        F.InterfaceReaction(
            id=1,
            subdomains=[vol_m, vol_l],
            k_plus=k_plus_HT,
            k_minus=k_minus,
            reactants=[H, T],
            products=[HT],
        ),
        F.InterfaceReaction(
            id=1,
            subdomains=[vol_m, vol_l],
            k_plus=k_plus_homo,
            k_minus=k_minus,
            reactants=[T, T],
            products=[T2],
        ),
    ]

    if duplicate_HT:
        model.interfaces.append(
            F.InterfaceReaction(
                id=1,
                subdomains=[vol_m, vol_l],
                k_plus=k_plus_HT,
                k_minus=k_minus,
                reactants=[T, H],
                products=[HT],
            )
        )

    # metal side held fixed; the liquid has no outlet, so the steady state is
    # a zero-flux state, i.e. equilibrium
    model.boundary_conditions = [
        F.FixedConcentrationBC(species=H, subdomain=left, value=c_H),
        F.FixedConcentrationBC(species=T, subdomain=left, value=c_T),
    ]

    model.temperature = 300
    model.settings = F.Settings(final_time=200, atol=1e-12, rtol=1e-12, stepsize=2.0)
    model.exports = [
        F.Profile1DExport(field=H2, subdomain=vol_l),
        F.Profile1DExport(field=HT, subdomain=vol_l),
        F.Profile1DExport(field=T2, subdomain=vol_l),
    ]

    model.initialise()
    model.run()

    # liquid is uniform, so take the interface node of each profile at each
    # exported time
    t = np.array(model.exports[0].t)
    series = np.array([[data[0] for data in e.data] for e in model.exports])
    return t, series


def plot(results, filename):
    """The three carriers relaxing to equilibrium, labelled on the curves.

    Left: equal forward constants, which give H2 : HT : T2 = 1 : 1 : 1. Right:
    the statistical factor k_HT = 2 k_homo, which gives 1 : 2 : 1. The
    duplicated-channel run is overlaid on the left panel to show that it only
    rescales the approach in time and lands on the same equilibrium.
    """
    mt.set_theme("urban")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    # the whole transient is over well before the end of the run; show only the
    # part where anything happens
    t_max = 30.0

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.2), sharey=True)

    for ax, key, case, ratio in zip(
        axes,
        ("plain", "statistical"),
        ("$k^+_{HT} = k^+_{HH}$", "$k^+_{HT} = 2\\,k^+_{HH}$"),
        ("1 : 1 : 1", "1 : 2 : 1"),
    ):
        t, series = results[key]
        window = t <= t_max
        t = t[window]
        c_H2, c_HT, c_T2 = series[:, window]

        # H2 and T2 coincide by symmetry, so draw H2 thick and pale underneath
        ax.plot(t, c_H2, color="C0", linewidth=3, alpha=0.4)
        ax.plot(t, c_T2, color="C1", linestyle="--")
        ax.plot(t, c_HT, color="C2")

        ax.annotate(
            "HT",
            xy=(t[-1], c_HT[-1]),
            xytext=(-22, 5),
            textcoords="offset points",
            color="C2",
            weight="bold",
        )
        ax.annotate(
            "H$_2$, T$_2$",
            xy=(t[-1], c_H2[-1]),
            xytext=(-58, -22),
            textcoords="offset points",
            color="C0",
            weight="bold",
        )
        ax.annotate(
            f"{case}\nH$_2$ : HT : T$_2$ = {ratio}",
            xy=(0.04, 0.80),
            xycoords="axes fraction",
            color="C0",
            weight="bold",
        )
        ax.set_xlabel("time")
        ax.set_xlim(0, t_max)

    # the duplication trick, on the panel whose equilibrium it fails to change
    t_dup, series_dup = results["duplicated"]
    window = t_dup <= t_max
    axes[0].plot(
        t_dup[window], series_dup[1][window], color="C2", linestyle=":", linewidth=1.5
    )
    axes[0].annotate(
        "HT, channel declared twice:\nfaster, same equilibrium",
        xy=(t_dup[3], series_dup[1][3]),
        xytext=(18, -62),
        textcoords="offset points",
        color="C2",
        weight="bold",
    )

    axes[0].set_ylabel("concentration")
    axes[0].set_ylim(0, 2.3)
    fig.tight_layout()
    fig.savefig(filename)


if __name__ == "__main__":
    cases = [
        ("plain", "plain mass action     (k_HT = k_homo)", 1.0, False),
        ("statistical", "statistical factor    (k_HT = 2 k_homo)", 2.0, False),
        ("duplicated", "duplicated channel    (H+T and T+H, k_HT = k_homo)", 1.0, True),
    ]
    results = {}
    for key, label, k_HT, duplicate in cases:
        t, series = run(k_plus_HT=k_HT, duplicate_HT=duplicate)
        results[key] = (t, series)
        c_H2, c_HT, c_T2 = series[:, -1]
        K = c_HT**2 / (c_H2 * c_T2)
        print(f"\n{label}")
        print(f"  c_H2 = {c_H2:.6f}  c_HT = {c_HT:.6f}  c_T2 = {c_T2:.6f}")
        print(f"  ratio H2 : HT : T2 = 1 : {c_HT / c_H2:.4f} : {c_T2 / c_H2:.4f}")
        print(f"  K_exch = {K:.6f}")

    plot(results, "isotopologue_equilibrium.pdf")
