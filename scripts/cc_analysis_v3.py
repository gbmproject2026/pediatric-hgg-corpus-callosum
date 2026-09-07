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


def wilson(k, n, z=1.96):
    if n == 0: return (0., 0.)
    p = k / n; den = 1 + z*z/n
    c = (p + z*z/(2*n))/den
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return (max(0., c-h)*100, min(1., c+h)*100)


def crude(a, b, c, dd):
    a_, b_, c_, d_ = a+.5, b+.5, c+.5, dd+.5
    l = math.log(a_*d_/(b_*c_)); se = math.sqrt(1/a_+1/b_+1/c_+1/d_)
    return math.exp(l), math.exp(l-1.96*se), math.exp(l+1.96*se)


def adj(pset, gset, col):
    f = pd.concat([pset.assign(g=1), gset.assign(g=0)])
    f = f[f.WT_vol > 0]
    m = sm.Logit(f[col], sm.add_constant(f[["g", "log_wt"]])).fit(disp=False)
    lo, hi = m.conf_int().loc["g"]
    return math.exp(m.params["g"]), math.exp(lo), math.exp(hi), float(m.pvalues["g"])


rows = []


def add(label, col, pset, primary=False):
    p_, g_ = pset.copy(), gbm.copy()
    p_[col] = p_[col].astype(int); g_[col] = g_[col].astype(int)
    a, na = int(p_[col].sum()), len(p_)
    b, nb = int(g_[col].sum()), len(g_)
    _, p = sps.fisher_exact([[a, na-a], [b, nb-b]])
    o, lo, hi = crude(a, na-a, b, nb-b)
    ao, alo, ahi, ap = adj(p_, g_, col)
    pl, ph = wilson(a, na); gl, gh = wilson(b, nb)
    rows.append(dict(label=label, primary=primary, pk=a, pn=na, ped=100*a/na,
                     ped_lo=pl, ped_hi=ph, gk=b, gn=nb, gbm=100*b/nb,
                     gbm_lo=gl, gbm_hi=gh, p=p, or_=o, lo=lo, hi=hi,
                     aor=ao, alo=alo, ahi=ahi, ap=ap))


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
print(f"{'outcome':34} {'pHGG':>16} {'GBM':>17} {'OR (95% CI)':>21} {'p':>7} {'p BH':>7} {'adjOR':>7}")
print("-" * 112)
for _, r in R.iterrows():
    star = " *" if r.primary else "  "
    print(f"{r.label:32}{star} {r.pk:>3}/{r.pn:<3} {r.ped:>5.1f}% {r.gk:>5}/{r.gn:<4} {r.gbm:>5.1f}% "
          f"{r.or_:>7.2f} ({r.lo:.2f}-{r.hi:>5.2f}) {ps(r.p):>7} {ps(r.p_bh):>7} {r.aor:>7.2f}")
print(f"\np BH = Benjamini-Hochberg across all {len(R)} outcomes. "
      f"Survives BH at .05: {int((R.p_bh < .05).sum())}/{len(R)}")

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
