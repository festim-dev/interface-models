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

import festim as F
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

    # liquid is uniform at equilibrium; take the interface node
    return tuple(e.data[-1][0] for e in model.exports)


if __name__ == "__main__":
    cases = [
        ("plain mass action     (k_HT = k_homo)", 1.0, False),
        ("statistical factor    (k_HT = 2 k_homo)", 2.0, False),
        ("duplicated channel    (H+T and T+H, k_HT = k_homo)", 1.0, True),
    ]
    for label, k_HT, duplicate in cases:
        c_H2, c_HT, c_T2 = run(k_plus_HT=k_HT, duplicate_HT=duplicate)
        K = c_HT**2 / (c_H2 * c_T2)
        print(f"\n{label}")
        print(f"  c_H2 = {c_H2:.6f}  c_HT = {c_HT:.6f}  c_T2 = {c_T2:.6f}")
        print(f"  ratio H2 : HT : T2 = 1 : {c_HT / c_H2:.4f} : {c_T2 / c_H2:.4f}")
        print(f"  K_exch = {K:.6f}")
