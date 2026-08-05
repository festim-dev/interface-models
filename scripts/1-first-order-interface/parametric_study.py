from example_usage import run_model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import morethemes as mt


def damkohler_number(D, k, L):
    return k * L / D


D_A = 0.5
D_B = 1.0

alpha = 0.5
k_plus = 3
k_minus = k_plus * alpha
x_interface = 0.5

c_0 = 2
c_L = 1


all_models = []
all_models_steady = []
all_damkohler_numbers = np.logspace(-2, 2, 10)
for da in all_damkohler_numbers:
    k_plus = da * D_A / x_interface
    my_model = run_model(D_A, D_B, k_plus, k_plus * alpha, x_interface, c_0, c_L)
    my_model_steady = run_model(
        D_A, D_B, k_plus, k_plus * alpha, x_interface, c_0, c_L, transient=False
    )
    all_models_steady.append(my_model_steady)
    all_models.append(my_model)

mt.set_theme("urban")

# https://coolors.co/1a4848-f7b000-f46036-c9f2c7-aceca1
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    color=["#1a4848", "#f7b000", "#f46036", "#c9f2c7", "#aceca1"]
)

# remove top and right spines for all plots in matplotlib params
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.figure()

# red to blue colormap
cmap = plt.get_cmap("RdYlGn")

# center colormap around 1
norm = LogNorm(
    vmin=min(all_damkohler_numbers), vmax=max(all_damkohler_numbers), clip=False
)
norm.vmin = min(all_damkohler_numbers)
norm.vmax = max(all_damkohler_numbers)

for model, da in zip(all_models, all_damkohler_numbers):
    all_ts = np.array(model.exports[0].t)
    all_ratios = []
    for idx, data in enumerate(model.exports[0].data):
        c_A_int = model.exports[0].data[idx][-1]
        c_B_int = model.exports[1].data[idx][0]
        ratio = c_A_int / c_B_int
        all_ratios.append(ratio)

    plt.plot(
        all_ts,
        all_ratios,
        color=cmap(norm(da)),
    )
plt.axhline(y=alpha, color="C0", linestyle="--")
plt.annotate(
    f"$k_-/k_+ = {alpha}$",
    xy=(0.5, alpha),
    xytext=(0.3, alpha - 0.08),
    color="C0",
    weight="bold",
)

plt.annotate("Da $\ll 1$", xy=(0.1, 1.2), color=cmap(norm(0.01)), weight="bold")
plt.annotate("Da $\gg 1$", xy=(0.45, alpha + 0.1), color=cmap(norm(100)), weight="bold")


plt.xlabel("time")
plt.ylabel("$c_A/c_B$")
# plt.ylim(bottom=0)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Required for ScalarMappable
cbar = plt.colorbar(
    sm,
    label="Damköhler Number (Da)",
    ax=plt.gca(),
)

plt.tight_layout()
plt.savefig("parametric_study_damkohler.pdf")
plt.show()

plt.figure(figsize=(6, 3))

# plot c_A/c_B as a function of Da
all_ratios = []
for model, da in zip(all_models_steady, all_damkohler_numbers):
    c_A_int = model.exports[0].data[-1][-1]
    c_B_int = model.exports[1].data[-1][0]
    ratio = c_A_int / c_B_int
    all_ratios.append(ratio)

plt.plot(
    all_damkohler_numbers,
    all_ratios,
    color=cmap(norm(da)),
    marker="o",
    markersize=4,
    alpha=0.5,
)

plt.axhline(y=alpha, color="C0", linestyle="--")
plt.annotate(
    "$k_-/k_+$",
    xy=(1, alpha),
    xytext=(0, -15),
    color="C0",
    textcoords="offset points",
    weight="bold",
)

plt.xlabel("Damköhler Number (Da)")
plt.ylabel("$c_A/c_B$")
plt.xscale("log")
plt.ylim(bottom=0)
plt.tight_layout()
plt.savefig("parametric_study_damkohler_steady.pdf")
plt.show()
