"""Hydrogen diffusivity in tungsten from the HTM compiled database (database.json release asset)."""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
db, out = sys.argv[1], sys.argv[2]
k_B = 8.617333e-5
props = [p for p in json.load(open(db))["data"] if p["material"] == "tungsten" and p["type"] == "Diffusivity"]
plt.rcParams.update({"font.family": "Inter", "font.size": 20, "axes.edgecolor": "#1a4848", "axes.labelcolor": "#1a4848",
                     "xtick.color": "#1a4848", "ytick.color": "#1a4848"})
fig, ax = plt.subplots(figsize=(8, 6), facecolor="#efefef")
ax.set_facecolor("#efefef")
for s in ["top", "right"]: ax.spines[s].set_visible(False)
cols = {"H": "#1a4848", "D": "#f7b000", "T": "#e4572e"}
vals_500 = []
for p in props:
    iso = (p.get("isotope") or "H").upper()
    D0, E = p["pre_exp"]["value"], p["act_energy"]["value"]
    Tmin, Tmax = p["range"]["value"] if p.get("range") else (300, 1200)
    T = np.linspace(Tmin, Tmax, 50)
    ax.semilogy(1000 / T, D0 * np.exp(-E / k_B / T), color=cols.get(iso, "#1a4848"), lw=2.2, alpha=0.85)
    vals_500.append(D0 * np.exp(-E / k_B / 500))
for iso, c in cols.items():
    ax.plot([], [], color=c, lw=3, label=iso)
ax.legend(frameon=False, loc="upper right", fontsize=18, ncol=3)
ax.set_xlabel("1000 / T (K$^{-1}$)")
ax.set_ylabel("diffusivity (m$^2$ s$^{-1}$)")
sec = ax.secondary_xaxis("top", functions=(lambda x: 1000 / np.maximum(x, 1e-9), lambda T: 1000 / np.maximum(T, 1e-9)))
sec.set_xticks([1000, 600, 400, 300, 250]); sec.set_xlabel("T (K)")
fig.tight_layout()
fig.savefig(out, dpi=150, facecolor="#efefef")
v = np.array(vals_500)
print(len(props), "records")
