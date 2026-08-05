from mpi4py import MPI

import numpy as np
import pyvista
from dolfinx import plot
from dolfinx.mesh import create_unit_square

import festim as F

import matplotlib.pyplot as plt
from analytical_solution import solution

import morethemes as mt


def run_model(D_A, D_B, k_plus, k_minus, x_interface, c_0, c_L, transient=True):
    assert 0 < x_interface < 1, "x_interface must be between 0 and 1"
    my_model = F.HydrogenTransportProblemDiscontinuous()

    mat1 = F.Material(D_0=D_A, E_D=0)
    mat2 = F.Material(D_0=D_B, E_D=0)

    vol1 = F.VolumeSubdomain1D(id=1, material=mat1, borders=[0, x_interface])
    vol2 = F.VolumeSubdomain1D(id=2, material=mat2, borders=[x_interface, 1])

    left = F.SurfaceSubdomain(id=3, locator=lambda x: np.isclose(x[0], 0))
    right = F.SurfaceSubdomain(id=4, locator=lambda x: np.isclose(x[0], 1))

    my_model.subdomains = [vol1, vol2, left, right]

    # dolfinx_mesh = create_unit_square(MPI.COMM_WORLD, 10, 10)
    # my_model.mesh = F.Mesh(dolfinx_mesh)

    my_model.mesh = F.Mesh1D(vertices=np.linspace(0, 1, 101))

    A = F.Species("A", subdomains=[vol1, vol2])

    my_model.species = [A]

    my_model.interfaces = [
        F.InterfaceReaction(
            id=1,
            subdomains=[vol1, vol2],
            k_plus=k_plus,
            k_minus=k_minus,
            reactants=[A],
            products=[A],
        ),
    ]

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(species=A, subdomain=left, value=c_0),
        F.FixedConcentrationBC(species=A, subdomain=right, value=c_L),
    ]

    my_model.temperature = 300

    my_model.settings = F.Settings(atol=1e-9, rtol=1e-9)
    if transient:
        my_model.settings.stepsize = 0.02
        my_model.settings.final_time = 0.5
    else:
        my_model.settings.transient = False

    my_model.exports = [
        F.Profile1DExport(field=A, subdomain=vol1),
        F.Profile1DExport(field=A, subdomain=vol2),
    ]

    from dolfinx.log import set_log_level, LogLevel

    set_log_level(LogLevel.INFO)

    my_model.initialise()
    my_model.run()
    return my_model


if __name__ == "__main__":
    D_A = 2
    D_B = 1

    k_plus = 3
    k_minus = 1
    x_interface = 0.6

    c_0 = 2
    c_L = 1

    my_model = run_model(D_A, D_B, k_plus, k_minus, x_interface, c_0, c_L)

    mt.set_theme("urban")

    # https://coolors.co/1a4848-f7b000-f46036-c9f2c7-aceca1
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
    )

    # remove top and right spines for all plots in matplotlib params
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False

    plt.figure(figsize=(6, 3))
    for e in my_model.exports:
        for idx, data in enumerate(e.data):
            plt.plot(
                e.x,
                e.data[idx],
                label=f"Subdomain {e.subdomain.id}",
                color="C0",
                alpha=0.5,
            )

    c_a_int, c_b_int, phi = solution(
        c_0=c_0,
        c_L=c_L,
        D_A=D_A,
        D_B=D_B,
        L_A=x_interface,
        L_B=1 - x_interface,
        k_plus=k_plus,
        k_minus=k_minus,
    )

    plt.plot(
        [0, x_interface, None, x_interface, 1],
        [c_0, c_a_int, None, c_b_int, c_L],
        "--",
        color="C1",
        label="Analytical solution",
    )
    plt.annotate(
        "Analytical solution",
        weight="bold",
        xy=(0.5, c_a_int),
        xytext=(0.6, 0.5),
        color="C1",
    )

    plt.ylim(bottom=0)
    plt.xlabel("x [m]")
    plt.ylabel("c [m$^{-3}$]")
    plt.tight_layout()
    plt.show()

    plt.figure()
    all_ts = np.array(my_model.exports[0].t)
    all_ratios = []
    for idx, data in enumerate(my_model.exports[0].data):
        c_A_int = my_model.exports[0].data[idx][-1]
        c_B_int = my_model.exports[1].data[idx][0]
        ratio = c_A_int / c_B_int
        all_ratios.append(ratio)

    plt.plot(all_ts, all_ratios, label="Numerical solution", color="C0")
    plt.axhline(y=k_minus / k_plus, color="C1", linestyle="--")
    plt.xlabel("t")
    plt.ylabel("$c_A/c_B$")
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()
