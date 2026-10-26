import csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
base = "/root/iaea/ndw/section_3_model_parameter_evaluation/data/"
area = 12e-3 * 15e-3
t_start = 72 * 3600 + 0.75 * 0.5 * 24 * 3600
dpas = [2.5, 0.5, 0.23, 0.1, 0.023, 0.005, 0.001, 0]
cmap = LinearSegmentedColormap.from_list("ngt", ["#f7b000", "#e4572e", "#1a4848"])
norm = LogNorm(1e-3, 2.5)
plt.rcParams.update({"font.size": 23, "axes.edgecolor": "#1a4848", "axes.labelcolor": "#1a4848",
                     "xtick.color": "#1a4848", "ytick.color": "#1a4848"})
fig, ax = plt.subplots(figsize=(9, 5.4), facecolor="#efefef")
ax.set_facecolor("#efefef")
for s in ["top", "right"]: ax.spines[s].set_visible(False)
for i, d in enumerate(dpas):
    exp = np.genfromtxt(f"{base}tds_data_schwartz_selinger/{d}_dpa.csv", delimiter=",")
    col = "#7a8a8a" if d == 0 else cmap(norm(d))
    ax.plot(exp[:, 0], exp[:, 1] / area / 1e17, color=col, lw=3)
    T, f1, f2 = [], [], []
    with open(f"{base}damaged_sample_tds_fittings/dpa_{d}/last.csv") as fh:
        for row in csv.reader(fh):
            if "t(s)" in row: continue
            if float(row[0]) >= t_start:
                T.append(float(row[1])); f1.append(float(row[2])); f2.append(float(row[3]))
    ax.plot(T, (-np.array(f1) - np.array(f2)) / 1e17, color="#1a4848", lw=1.6, ls="--", alpha=0.8,
            label="FESTIM" if i == 0 else None)
ax.plot([], [], color="#e4572e", lw=3, label="experiment")
ax.set_xlim(300, 1000); ax.set_ylim(0, 1)
ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"flux ($10^{17}$ D m$^{-2}$ s$^{-1}$)")
ax.legend(frameon=False, loc="upper right")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cb = fig.colorbar(sm, ax=ax, pad=0.02); cb.set_label("damage (dpa)"); cb.outline.set_visible(False)
ax.plot([], [], color="#7a8a8a", lw=3, label="undamaged")
ax.legend(frameon=False, loc="upper right", fontsize=20)
fig.tight_layout()
fig.savefig("/root/iaea/web/neutron_damage_tds.png", dpi=150, facecolor="#efefef")
