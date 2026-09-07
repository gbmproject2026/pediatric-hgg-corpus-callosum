"""cc_figures_v3.py  Figures on the corrected measurements (JHU atlas, fixed labels)."""
import math
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import config as cfg
OUT = cfg.FIGURES
OUT.mkdir(parents=True, exist_ok=True)
R = pd.read_csv(cfg.OUTCOMES)
PED, GBM_, GREY = "#c1121f", "#2d6a9f", "#9aa0a6"
def ps(p): return "p<.001" if p < .001 else f"p={p:.3f}".replace("0.", ".")

# FIG A: forest, unadjusted primary
fig = plt.figure(figsize=(13.8, 8.2))
gs = fig.add_gridspec(1, 2, width_ratios=[3.1, 1.0], wspace=.30)
ax = fig.add_subplot(gs[0, 0])
y = np.arange(len(R))[::-1].astype(float)
for i, r in R.iterrows():
    sig = r.p < .05
    col = PED if (sig and r.pr > 1) else (GBM_ if (sig and r.pr < 1) else GREY)
    ax.plot([r.pr_lo, r.pr_hi], [y[i]]*2, color=col, lw=2.7, solid_capstyle="round")
    ax.plot(r.pr, y[i], "o", color=col, ms=8.5, zorder=3,
            mec="black" if r.primary else col, mew=1.7 if r.primary else 0)
ax.axvline(1, color="black", ls=":", lw=1.5)
ax.set_yticks(y); ax.set_yticklabels([r.label + ("  ★" if r.primary else "") for _, r in R.iterrows()], fontsize=9.5)
ax.set_xscale("log"); ax.set_xlim(.3, 12)
ax.set_xlabel("Prevalence ratio, pediatric HGG vs adult GBM (log scale)", fontsize=10.5)
for i, r in R.iterrows():
    ax.text(11.5, y[i], f"{r.ped:.1f}% vs {r.gbm:.1f}%   {ps(r.p)}", fontsize=8.2, va="center", ha="right",
            color="#222" if r.p < .05 else "#8a8a8a")
ax.grid(alpha=.22, axis="x")
ax.text(.33, -0.95, "more common in adult GBM", fontsize=8.5, color=GBM_, style="italic")
ax.text(2.6, -0.95, "more common in pediatric HGG", fontsize=8.5, color=PED, style="italic")
ax.set_title("Corpus callosum involvement and midline crossing\n"
             "JHU ICBM-DTI-81 atlas; pHGG n=114 (80 enhancing) vs adult GBM n=1213", fontsize=12)
ax2 = fig.add_subplot(gs[0, 1])
sel = R[R.label.isin(["CC involvement, TC", "Midline crossing, TC (>20/side)"])].reset_index(drop=True)
yy = [1., 0.]
for k, r in sel.iterrows():
    ax2.plot([r.pr_lo, r.pr_hi], [yy[k]+.13]*2, color=GREY, lw=2.6, solid_capstyle="round")
    ax2.plot(r.pr, yy[k]+.13, "o", color=GREY, ms=8)
    ax2.plot([r.pr_adj_lo, r.pr_adj_hi], [yy[k]-.13]*2, color="#7b2d26", lw=2.6, solid_capstyle="round")
    ax2.plot(r.pr_adj, yy[k]-.13, "s", color="#7b2d26", ms=8)
    ax2.text(r.pr_hi*1.1, yy[k]+.13, f"{r.pr:.2f}", fontsize=8, va="center", color="#666")
    ax2.text(r.pr_adj_hi*1.1, yy[k]-.13, f"{r.pr_adj:.2f}", fontsize=8, va="center", color="#7b2d26")
ax2.axvline(1, color="black", ls=":", lw=1.4)
ax2.set_yticks(yy); ax2.set_yticklabels(["CC involvement\ntumor core", "Midline crossing\ntumor core"], fontsize=9)
ax2.set_xscale("log"); ax2.set_xlim(.8, 8); ax2.set_ylim(-.6, 1.6)
from matplotlib.ticker import FixedLocator, NullLocator, FixedFormatter
ax2.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8]))
ax2.xaxis.set_major_formatter(FixedFormatter(["1", "2", "4", "8"]))
ax2.xaxis.set_minor_locator(NullLocator())
ax2.set_xlabel("Prevalence ratio", fontsize=9.5); ax2.grid(alpha=.22, axis="x")
ax2.set_title("Sensitivity analysis\nvolume adjustment", fontsize=10.5)
ax2.legend(handles=[Line2D([], [], color=GREY, marker="o", lw=2.4, label="unadjusted"),
                    Line2D([], [], color="#7b2d26", marker="s", lw=2.4, label="volume-adjusted")],
           loc="upper right", fontsize=8, framealpha=.95)
fig.text(.685, .045, "Adult GBM tumors were 2.45x larger (p<.001). Volume plausibly mediates the\n"
                     "tumor-type/callosal relationship, so unadjusted estimates are primary.\n"
                     "Adjusted values are modified Poisson prevalence ratios.",
         fontsize=7.4, color="#555", va="top")
fig.savefig(OUT/"Figure4_forest.png", dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)

# FIG B: subregion profile, WT and TC
d = pd.read_csv(cfg.SIDES)
ped, gbm = d[d.cohort=="pHGG"], d[d.cohort=="GBM"]
pe = ped[ped.ET_vol>0]
fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
for ax, comp, pset, ttl in [(axes[0], "WT", ped, "Whole tumor (all cases)"),
                            (axes[1], "TC", pe, "Tumor core (enhancing cases)")]:
    regs = ["genu", "body", "splenium"]
    pv = [100*(pset[f"{r}_{comp}_n"]>0).mean() for r in regs]
    gv = [100*(gbm[f"{r}_{comp}_n"]>0).mean() for r in regs]
    x = np.arange(3); w = .35
    ax.bar(x-w/2, pv, w, color=PED, label=f"pediatric HGG (n={len(pset)})", zorder=3)
    ax.bar(x+w/2, gv, w, color=GBM_, label=f"adult GBM (n={len(gbm)})", zorder=3)
    import scipy.stats as sps
    for i, r in enumerate(regs):
        a=int((pset[f"{r}_{comp}_n"]>0).sum()); b=int((gbm[f"{r}_{comp}_n"]>0).sum())
        _, p = sps.fisher_exact([[a,len(pset)-a],[b,len(gbm)-b]])
        ax.text(x[i], max(pv[i],gv[i])+3.5, ps(p), ha="center", fontsize=8.5, style="italic")
        ax.text(x[i]-w/2, pv[i]+.8, f"{pv[i]:.1f}%", ha="center", fontsize=8.5, fontweight="bold")
        ax.text(x[i]+w/2, gv[i]+.8, f"{gv[i]:.1f}%", ha="center", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(["Genu","Body","Splenium"], fontsize=11)
    ax.set_ylabel("% of cohort involved", fontsize=10.5); ax.set_ylim(0, 92)
    ax.set_title(ttl, fontsize=11); ax.legend(fontsize=8.5); ax.grid(alpha=.22, axis="y")
fig.suptitle("Subregional corpus callosum involvement, anterior to posterior", fontsize=13, y=1.0)
fig.tight_layout(); fig.savefig(OUT/"Figure3_subregions.png", dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)

# FIG C: prevalence by compartment
sel = R[R.label.str.contains("CC involvement|Midline crossing")].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(11.5, 5.8))
x = np.arange(len(sel)); w=.35
ax.bar(x-w/2, sel.ped, w, color=PED, label="pediatric HGG", zorder=3)
ax.bar(x+w/2, sel.gbm, w, color=GBM_, label="adult GBM", zorder=3)
for i, r in sel.iterrows():
    ax.text(x[i]-w/2, r.ped+1.5, f"{r.ped:.1f}%", ha="center", fontsize=8.3, fontweight="bold")
    ax.text(x[i]+w/2, r.gbm+1.5, f"{r.gbm:.1f}%", ha="center", fontsize=8.3)
    ax.text(x[i], max(r.ped,r.gbm)+9, ps(r.p), ha="center", fontsize=8, style="italic")
    ax.text(x[i], -11, f"{r.pk}/{r.pn}\n{r.gk}/{r.gn}", ha="center", fontsize=7, color="#555")
ax.set_xticks(x); ax.set_xticklabels([s.replace(", ","\n").replace(" (>20/side)","").replace(" (>5/side)","") for s in sel.label], fontsize=8.6)
ax.set_ylabel("% of cohort", fontsize=10.5); ax.set_ylim(-16, 100); ax.axhline(0, color="black", lw=.8)
ax.set_title("Corpus callosum involvement and midline crossing by compartment\n"
             "counts below axis (pHGG / GBM)", fontsize=11.5)
ax.legend(fontsize=9); ax.grid(alpha=.22, axis="y")
fig.tight_layout(); fig.savefig(OUT/"Figure2_prevalence.png", dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
print("wrote Figure2_prevalence.png, Figure3_subregions.png, Figure4_forest.png ->", OUT)
