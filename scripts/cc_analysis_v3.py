"""
cc_analysis_v3.py
Definitive analysis on cc_sides_v3.csv (JHU atlas, corrected labels).

Primary estimates are UNADJUSTED (tumor volume is a mediator, not a confounder).
Volume-adjusted logistic regression is reported as a sensitivity analysis.
"""
import math
from pathlib import Path
import numpy as np, pandas as pd, scipy.stats as sps, statsmodels.api as sm
import config as cfg
d = pd.read_csv(cfg.SIDES)
ped, gbm = d[d.cohort == "pHGG"].copy(), d[d.cohort == "GBM"].copy()
ped_enh = ped[ped.ET_vol > 0].copy()
for x in (ped, gbm, ped_enh):
    x.loc[:, "log_wt"] = np.log10(x.WT_vol.clip(lower=1))

CROSS_THR = {"WT": 20, "TC": 20, "ET": 5}


def prev_ratio(k1, n1, k0, n0):
    """Prevalence ratio with 95% CI by the Katz log method.

    The design is cross-sectional and the outcomes are common (18% to 84%), so a
    prevalence ratio is directly interpretable. An odds ratio would overstate the
    difference substantially at these prevalences.
    """
    a, b = k1 + .5, k0 + .5           # continuity correction for empty cells
    p1, p0 = a / (n1 + 1), b / (n0 + 1)
    pr = p1 / p0
    se = math.sqrt((1 - p1) / a + (1 - p0) / b)
    return pr, pr * math.exp(-1.96 * se), pr * math.exp(1.96 * se)


def adj(pset, gset, col):
    """Volume-adjusted prevalence ratio by modified Poisson regression.

    Poisson regression with robust (HC0) standard errors, which is the standard
    way to obtain an adjusted prevalence ratio for a common binary outcome.
    Logistic regression would return an odds ratio and overstate the effect.
    """
    f = pd.concat([pset.assign(g=1), gset.assign(g=0)])
    f = f[f.WT_vol > 0]
    m = sm.GLM(f[col].astype(int), sm.add_constant(f[["g", "log_wt"]]),
               family=sm.families.Poisson()).fit(cov_type="HC0")
    lo, hi = m.conf_int().loc["g"]
    return math.exp(m.params["g"]), math.exp(lo), math.exp(hi), float(m.pvalues["g"])


rows = []


def add(label, col, pset, primary=False):
    p_, g_ = pset.copy(), gbm.copy()
    p_[col] = p_[col].astype(int); g_[col] = g_[col].astype(int)
    a, na = int(p_[col].sum()), len(p_)
    b, nb = int(g_[col].sum()), len(g_)
    _, p = sps.fisher_exact([[a, na-a], [b, nb-b]])
    o, lo, hi = prev_ratio(a, na, b, nb)
    ao, alo, ahi, ap = adj(p_, g_, col)
    rows.append(dict(label=label, primary=primary, pk=a, pn=na, ped=100*a/na,
                     gk=b, gn=nb, gbm=100*b/nb, p=p,
                     pr=o, pr_lo=lo, pr_hi=hi,
                     pr_adj=ao, pr_adj_lo=alo, pr_adj_hi=ahi, p_adj=ap))


for r in ["whole", "genu", "body", "splenium"]:
    for c in ["WT", "TC", "ET"]:
        col = f"inv_{r}_{c}"
        for x in (ped, gbm, ped_enh):
            x.loc[:, col] = x[f"{r}_{c}_n"] > 0
        pset = ped_enh if c in ("TC", "ET") else ped
        nm = {"whole": "CC", "genu": "Genu", "body": "Body", "splenium": "Splenium"}[r]
        add(f"{nm} involvement, {c}", col, pset,
            primary=(r == "whole" and c == "TC"))

for c in ["WT", "TC", "ET"]:
    t = CROSS_THR[c]
    col = f"cross_{c}"
    for x in (ped, gbm, ped_enh):
        x.loc[:, col] = (x[f"whole_{c}_L"] > t) & (x[f"whole_{c}_R"] > t)
    pset = ped_enh if c in ("TC", "ET") else ped
    add(f"Midline crossing, {c} (>{t}/side)", col, pset)

R = pd.DataFrame(rows)


def bh_fdr(pv):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pv, float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


# Tumor core midline crossing is the pre-specified primary outcome and is
# reported unadjusted. The remaining outcomes are secondary, so a
# Benjamini-Hochberg adjusted p-value is reported across the family of 15.
R["p_bh"] = bh_fdr(R.p.values)
R.to_csv(cfg.OUTCOMES, index=False)


def ps(p): return "<.001" if p < 0.001 else f"{p:.3f}".lstrip("0")


print("=" * 108)
print("CORRECTED ANALYSIS  (JHU corpus callosum atlas; pHGG ET=1 TC=1+2; GBM ET=3 TC=1+3)")
print(f"pHGG n={len(ped)} (enhancing {len(ped_enh)}, non-enhancing {len(ped)-len(ped_enh)})   GBM n={len(gbm)}")
print("TC and ET comparisons use the enhancing pHGG subgroup; WT uses the full cohort.")
print("=" * 108)
print(f"{'outcome':34} {'pHGG':>16} {'GBM':>17} {'PR (95% CI)':>21} {'p':>7} {'p BH':>7} {'adjPR':>7}")
print("-" * 112)
for _, r in R.iterrows():
    star = " *" if r.primary else "  "
    print(f"{r.label:32}{star} {r.pk:>3}/{r.pn:<3} {r.ped:>5.1f}% {r.gk:>5}/{r.gn:<4} {r.gbm:>5.1f}% "
          f"{r.pr:>7.2f} ({r.pr_lo:.2f}-{r.pr_hi:>5.2f}) {ps(r.p):>7} {ps(r.p_bh):>7} {r.pr_adj:>7.2f}")
print(f"\np BH = Benjamini-Hochberg across all {len(R)} outcomes. "
      f"Survives BH at .05: {int((R.p_bh < .05).sum())}/{len(R)}")

# Continuous outcome: fraction of the corpus callosum involved. Avoids the
# arbitrary "any voxel" threshold and retains magnitude information.
CC_VOX = 35291
print(f"\n{'='*108}\nCONTINUOUS: fraction of the corpus callosum involved (Mann-Whitney U)\n{'='*108}")
print(f"{'measure':26} {'pHGG median':>13} {'GBM median':>12} {'p':>10}   restricted to involved cases")
for lab, col, pset in [("whole tumor", "whole_WT_n", ped), ("tumor core", "whole_TC_n", ped_enh)]:
    a = pset[col] / CC_VOX * 100
    b = gbm[col] / CC_VOX * 100
    uu = sps.mannwhitneyu(a, b)
    ai, bi = a[a > 0], b[b > 0]
    ui = sps.mannwhitneyu(ai, bi)
    print(f"{lab:26} {a.median():12.1f}% {b.median():11.1f}% {ps(uu.pvalue):>10}   "
          f"{ai.median():.1f}% vs {bi.median():.1f}%, p={ps(ui.pvalue)}")

u = sps.mannwhitneyu(ped_enh.WT_vol, gbm.WT_vol)
print(f"\nmedian WT volume: pHGG-all {ped.WT_vol.median():.0f}, pHGG-enh {ped_enh.WT_vol.median():.0f}, "
      f"GBM {gbm.WT_vol.median():.0f} ({gbm.WT_vol.median()/ped_enh.WT_vol.median():.2f}x, MWU p={u.pvalue:.1e})")

non = ped[ped.ET_vol == 0]
print(f"\nNon-enhancing pHGG (n={len(non)}):")
for r in ["whole", "genu", "body", "splenium"]:
    k = int((non[f"{r}_WT_n"] > 0).sum())
    print(f"  {r:9} WT involvement {k}/{len(non)} ({100*k/len(non):.1f}%)")
k = int(((non.whole_WT_L > 20) & (non.whole_WT_R > 20)).sum())
print(f"  crossing  WT {k}/{len(non)} ({100*k/len(non):.1f}%)")

print("\nSubregion distribution among CC-involved cases (share of involved CC voxels):")
for nm, x in [("pHGG", ped), ("GBM", gbm)]:
    inv = x[x.whole_WT_n > 0]
    tot = inv[["genu_WT_n", "body_WT_n", "splenium_WT_n"]].sum()
    print(f"  {nm:5} genu {100*tot.genu_WT_n/tot.sum():4.1f}%  "
          f"body {100*tot.body_WT_n/tot.sum():4.1f}%  splenium {100*tot.splenium_WT_n/tot.sum():4.1f}%")
print("\nWrote cc_results_v3.csv")
