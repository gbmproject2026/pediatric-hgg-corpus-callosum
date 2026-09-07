"""cc_figS6_reproducibility.py  Supplemental Figure S6: outcome agreement across\nregistration seeds. Reports Cohen kappa and per-seed prevalence, the measures that\ncorrespond to the binary outcomes actually reported."""
from pathlib import Path
import numpy as np, pandas as pd, scipy.stats as sps
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as cfg
SUP = cfg.FIGURES_SUP
D = pd.read_csv(SUP/"TableS4_reproducibility_per_case.csv")
S = pd.read_csv(SUP/"TableS6_outcome_stability.csv")

fig = plt.figure(figsize=(15.5, 5.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.9, 1.0], wspace=.24)
COL = {"WT": "#2d6a9f", "TC": "#f4a261", "ET": "#c1121f"}

# D: prevalence per seed
ax = fig.add_subplot(gs[0, 0])
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
ax.set_title("Outcome prevalence across independent registrations\n"
             "maximum shift 3.5 percentage points", fontsize=10.5)
ax.legend(fontsize=8.4, ncol=4, loc="upper right"); ax.grid(alpha=.25, axis="y")

# E: kappa summary
ax = fig.add_subplot(gs[0, 1])
y = np.arange(len(S))[::-1]
ax.barh(y, S.mean_kappa, .6, color="#2a9d8f", alpha=.88)
for i, r in S.iterrows():
    ax.text(r.mean_kappa+.008, y[i], f"{r.mean_kappa:.3f}", va="center", fontsize=8.6)
ax.axvline(.8, color="crimson", ls=":", lw=1.4)
ax.text(.805, len(S)-.6, "0.80\nsubstantial", fontsize=7.6, color="crimson")
ax.set_yticks(y); ax.set_yticklabels(S.outcome, fontsize=8.4)
ax.set_xlim(0, 1.05); ax.set_xlabel("Cohen kappa vs reference")
ax.set_title("Classification agreement", fontsize=10.5); ax.grid(alpha=.25, axis="x")

fig.suptitle("Supplemental Figure S6. Outcome agreement across three independent registration seeds "
             f"(pediatric cohort, n={len(D)})", fontsize=12.5, y=.98)
p = SUP/"FigureS6_reproducibility.png"
fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
print("wrote", p)
