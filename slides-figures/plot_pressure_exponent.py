"""Pressure exponent of the permeating flux at a metal/salt interface with H2 and HF channels.
Data: FESTIM sweep and analytical prediction from festim-dev/interface-models
(scripts/3-competing-channels/verification_model3_exponent.csv).
With Sieverts on the metal side, c_m = K_S sqrt(P), so J ~ P^(n/2)."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv, out = sys.argv[1], sys.argv[2]
rows = np.loadtxt(csv, delimiter=",", skiprows=1)
plt.rcParams.update({"font.family": "Inter", "font.size": 20, "axes.edgecolor": "#1a4848",
                     "axes.labelcolor": "#1a4848", "xtick.color": "#1a4848", "ytick.color": "#1a4848"})
fig, ax = plt.subplots(figsize=(8.4, 5.6), facecolor="#efefef")
ax.set_facecolor("#efefef")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

for level in (1.0, 0.5):
    ax.axhline(level, color="#8a9696", ls="--", lw=1.4, zorder=0)
for y, label, va in ((1.0, r"H$_2$ wins:  $J \propto P$", "bottom"), (0.5, r"HF wins:  $J \propto \sqrt{P}$", "top")):
    ax.annotate(label, xy=(1.5e-6, y), xytext=(0, 8 if va == "bottom" else -8), textcoords="offset points",
                va=va, ha="left", color="#5f6d6d", fontsize=18, weight="bold")

colours = {0.1: "#1a4848", 1.0: "#f7b000", 10.0: "#e4572e"}
names = {0.1: "reducing", 1.0: "", 10.0: "oxidising"}
for a_F, col in colours.items():
    sel = rows[:, 0] == a_F
    P = rows[sel, 2] ** 2
    m_pred, m_meas = rows[sel, 10] / 2, rows[sel, 9] / 2
    ax.semilogx(P, m_pred, color=col, lw=3)
    ax.semilogx(P[::4], m_meas[::4], "o", ms=8, mfc="none", mew=2, color=col)
    k = np.argmin(np.abs(m_pred - 0.75))
    lab = f"$a_F$ = {a_F:g}" + (f"\n{names[a_F]}" if names[a_F] else "")
    ax.annotate(lab, xy=(P[k], m_pred[k]), xytext=(10, -6), textcoords="offset points",
                color=col, fontsize=18, weight="bold", va="top")

ax.plot([], [], color="#1a4848", lw=3, label="analytical")
ax.plot([], [], "o", ms=8, mfc="none", mew=2, color="#1a4848", label="FESTIM")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.88), fontsize=17)
ax.set_ylim(0.44, 1.06)
ax.set_yticks([0.5, 0.75, 1.0], ["0.5", "0.75", "1"])
ax.set_xlim(1e-6, 1.5e2)
ax.set_xlabel(r"H$_2$ pressure at the interface (arb.)")
ax.set_ylabel("pressure exponent of the flux")
fig.tight_layout()
fig.savefig(out, dpi=150, facecolor="#efefef")
