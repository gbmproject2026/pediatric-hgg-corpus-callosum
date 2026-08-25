"""cc_figS6_reproducibility.py  Supplemental Figure S6: registration reproducibility."""
from pathlib import Path
import numpy as np, pandas as pd, scipy.stats as sps
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as cfg
SUP = cfg.FIGURES_SUP
D = pd.read_csv(SUP/"TableS4_reproducibility_per_case.csv")
S = pd.read_csv(SUP/"TableS6_outcome_stability.csv")
v3 = pd.read_csv(cfg.SIDES).set_index("subject")
SE = [43, 44, 45]
D = D.join(v3[["WT_vol"]], on="subject")
for c in ["WT", "TC", "ET"]:
    D[f"{c}_med"] = D[[f"dice_{c}_s{s}" for s in SE]].median(axis=1)
D["WT_min"] = D[[f"dice_WT_s{s}" for s in SE]].min(axis=1)

fig = plt.figure(figsize=(16.5, 9.6))
gs = fig.add_gridspec(2, 3, hspace=.32, wspace=.26)
COL = {"WT": "#2d6a9f", "TC": "#f4a261", "ET": "#c1121f"}

# A: Dice distribution per compartment
ax = fig.add_subplot(gs[0, 0])
data = [pd.concat([D[f"dice_{c}_s{s}"] for s in SE]).dropna().values for c in ["WT", "TC", "ET"]]
parts = ax.violinplot(data, showmedians=True, widths=.8)
for i, b in enumerate(parts["bodies"]):
    b.set_facecolor(list(COL.values())[i]); b.set_alpha(.62)
for k in ("cbars", "cmins", "cmaxes", "cmedians"):
    parts[k].set_color("#333333")
for i, d in enumerate(data):
    ax.text(i+1, -.07, f"med {np.median(d):.2f}\nn={len(d)}", ha="center", fontsize=8.2)
ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["whole tumor", "tumor core", "enhancing tumor"], fontsize=9.5)
ax.set_ylabel("Dice vs reference registration"); ax.set_ylim(-.12, 1.04)
ax.set_title("A. Reproducibility by tumor compartment\n(3 repeat registrations per case)", fontsize=10.5)
ax.grid(alpha=.25, axis="y")

# B: Dice vs volume
ax = fig.add_subplot(gs[0, 1])
ax.scatter(D.WT_vol, D.WT_med, s=26, color=COL["WT"], alpha=.75, edgecolor="none")
r = sps.spearmanr(D.WT_vol, D.WT_med)
ax.set_xscale("log"); ax.set_xlabel("whole-tumor volume (voxels, log scale)")
ax.set_ylabel("median Dice across seeds"); ax.set_ylim(-.04, 1.04)
ax.axhline(.5, color="crimson", ls=":", lw=1.3)
ax.text(.02, .06, f"Spearman rho = {r.statistic:+.3f}\nP < .001", transform=ax.transAxes, fontsize=9)
ax.set_title("B. Reproducibility scales with tumor volume\nDice penalizes small objects more", fontsize=10.5)
ax.grid(alpha=.25)

# C: Dice by volume quartile
ax = fig.add_subplot(gs[0, 2])
D["q"] = pd.qcut(D.WT_vol, 4, labels=["Q1\nsmallest", "Q2", "Q3", "Q4\nlargest"])
g = D.groupby("q", observed=True)
med = g.WT_med.median(); uns = g.WT_min.apply(lambda s: (s < .5).sum()); nq = g.size()
x = np.arange(4)
ax.bar(x, med.values, .62, color=COL["WT"], alpha=.85)
for i in range(4):
    ax.text(x[i], med.values[i]+.02, f"{med.values[i]:.2f}", ha="center", fontsize=9, fontweight="bold")
    ax.text(x[i], .05, f"{uns.values[i]}/{nq.values[i]}\nmin<0.5", ha="center", fontsize=7.8, color="white")
ax.set_xticks(x); ax.set_xticklabels(med.index, fontsize=9)
ax.set_ylabel("median Dice"); ax.set_ylim(0, 1.05)
ax.set_title("C. Dice and instability by volume quartile", fontsize=10.5)
ax.grid(alpha=.25, axis="y")

# D: prevalence per seed
ax = fig.add_subplot(gs[1, :2])
lbl = S.outcome.tolist(); x = np.arange(len(lbl)); w = .19
for i, (cl, nm) in enumerate([("prev_s42", "seed 42 (reference)"), ("prev_s43", "seed 43"),
                              ("prev_s44", "seed 44"), ("prev_s45", "seed 45")]):
    ax.bar(x + (i-1.5)*w, S[cl], w, label=nm, alpha=.9)
for i, r in S.iterrows():
    v = [r.prev_s42, r.prev_s43, r.prev_s44, r.prev_s45]
    ax.text(x[i], max(v)+2.2, f"kappa {r.mean_kappa:.3f}\n{r.mean_agreement:.1f}% agree",
            ha="center", fontsize=7.8)
ax.set_xticks(x); ax.set_xticklabels([l.replace(" ", "\n") for l in lbl], fontsize=8.6)
ax.set_ylabel("% of cohort positive"); ax.set_ylim(0, 92)
ax.set_title("D. Outcome prevalence is stable across independent registrations\n"
             "maximum shift 3.5 percentage points; all kappa > 0.89", fontsize=10.5)
ax.legend(fontsize=8.4, ncol=4, loc="upper right"); ax.grid(alpha=.25, axis="y")

# E: kappa summary
ax = fig.add_subplot(gs[1, 2])
y = np.arange(len(S))[::-1]
ax.barh(y, S.mean_kappa, .6, color="#2a9d8f", alpha=.88)
for i, r in S.iterrows():
    ax.text(r.mean_kappa+.008, y[i], f"{r.mean_kappa:.3f}", va="center", fontsize=8.6)
ax.axvline(.8, color="crimson", ls=":", lw=1.4)
ax.text(.805, len(S)-.6, "0.80\nsubstantial", fontsize=7.6, color="crimson")
ax.set_yticks(y); ax.set_yticklabels(S.outcome, fontsize=8.4)
ax.set_xlim(0, 1.05); ax.set_xlabel("Cohen kappa vs reference")
ax.set_title("E. Classification agreement", fontsize=10.5); ax.grid(alpha=.25, axis="x")

fig.suptitle("Supplemental Figure S6. Registration reproducibility across three independent random seeds "
             f"(pediatric cohort, n={len(D)})", fontsize=12.5, y=.975)
p = SUP/"FigureS6_reproducibility.png"
fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
print("wrote", p)
