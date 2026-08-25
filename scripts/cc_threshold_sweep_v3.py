"""
cc_threshold_sweep_v3.py
Which midline-crossing definition is best?

Two families are compared on the JHU atlas with corrected labels:
  ABSOLUTE  crossing if > T voxels of the compartment in the CC on each side
  FRACTIONAL crossing if > F% of each hemisphere's CC volume is involved
The JHU CC mask is 7.9% asymmetric about the split plane (left 18,341 vs right
16,950 voxels), so an absolute threshold is systematically easier to clear on the
left. A fractional threshold removes that bias.

Definitions are judged on THREE criteria, deliberately not on p-value alone:
  1. reproducibility  Cohen kappa of the classification across three independent
                      registration seeds (pediatric cohort)
  2. stability        whether cohort discrimination sits on a plateau or a spike
  3. discrimination   Fisher p and odds ratio, reported but NOT used to select
"""
import math
from pathlib import Path
import ants, nibabel as nib, numpy as np, pandas as pd, scipy.stats as sps
import config as cfg
OUTD = cfg.FIGURES_SUP
CC = ants.image_read(str(cfg.ATLAS/"CC_mask.nii.gz")).numpy() > 0
MID = CC.shape[0]//2
NL, NR = int(CC[:MID].sum()), int(CC[MID:].sum())
COH = open(cfg.COHORT).read().split()
SEEDS = [43, 44, 45]
LAB = {"pHGG": {"WT": None, "TC": (1, 2, 3), "ET": (1,)},
       "GBM":  {"WT": None, "TC": (1, 3), "ET": (3,)}}

def sides(seg, cohort, comp):
    labs = LAB[cohort][comp]
    arr = (seg > 0) if labs is None else np.isin(seg, labs)
    c = arr & CC
    return int(c[:MID].sum()), int(c[MID:].sum())

def collect(cohort, ids, regdir):
    out = {}
    for sid in ids:
        p = Path(regdir)/f"{sid}_seg_MNI.nii.gz"
        if not p.exists(): continue
        seg = np.asanyarray(nib.load(p).dataobj)
        out[sid] = {c: sides(seg, cohort, c) for c in ["WT", "TC", "ET"]}
    return out

import glob, os
gli_ids = [os.path.basename(f).replace("_seg_MNI.nii.gz", "")
           for f in glob.glob(str(cfg.GLI_REG/"*_seg_MNI.nii.gz"))]
ref_p = collect("pHGG", COH, cfg.PEDS_REG)
ref_g = collect("GBM", gli_ids, cfg.GLI_REG)
seed_p = {s: collect("pHGG", COH, cfg.seed_reg(s)) for s in SEEDS}
v3 = pd.read_csv(cfg.SIDES).set_index("subject")
enh = set(v3[(v3.cohort == "pHGG") & (v3.ET_vol > 0)].index)
print(f"loaded pHGG {len(ref_p)}, GBM {len(ref_g)}, seeds {[len(seed_p[s]) for s in SEEDS]}")

def cross_abs(lr, t): return lr[0] > t and lr[1] > t
def cross_frac(lr, f): return (lr[0]/NL)*100 > f and (lr[1]/NR)*100 > f

def kappa(a, b):
    a, b = np.array(a, int), np.array(b, int)
    po = (a == b).mean()
    pe = a.mean()*b.mean() + (1-a.mean())*(1-b.mean())
    return (po-pe)/(1-pe) if pe < 1 else 1.0

rows = []
GRID = {"absolute": [0, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200],
        "fractional": [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]}
for comp in ["WT", "TC", "ET"]:
    pids = [s for s in ref_p if (s in enh if comp in ("TC", "ET") else True)]
    for fam, grid in GRID.items():
        fn = cross_abs if fam == "absolute" else cross_frac
        for t in grid:
            fp = [fn(ref_p[s][comp], t) for s in pids]
            fg = [fn(ref_g[s][comp], t) for s in ref_g]
            a, na, b, nb = sum(fp), len(fp), sum(fg), len(fg)
            _, p = sps.fisher_exact([[a, na-a], [b, nb-b]])
            A, B, C, D = a+.5, na-a+.5, b+.5, nb-b+.5
            o = math.exp(math.log(A*D/(B*C)))
            ks = []
            for s in SEEDS:
                common = [x for x in pids if x in seed_p[s]]
                ks.append(kappa([fn(ref_p[x][comp], t) for x in common],
                                [fn(seed_p[s][x][comp], t) for x in common]))
            rows.append(dict(compartment=comp, family=fam, threshold=t,
                             ped_pct=100*a/na, gbm_pct=100*b/nb, n_ped=na,
                             odds_ratio=o, p=p, kappa=float(np.mean(ks))))
S = pd.DataFrame(rows)
S.to_csv(OUTD/"TableS8_crossing_definition_sweep.csv", index=False)

def ps(p): return "<.001" if p < .001 else f"{p:.3f}".lstrip("0")
for comp in ["WT", "TC", "ET"]:
    for fam in ["absolute", "fractional"]:
        s = S[(S.compartment == comp) & (S.family == fam)]
        unit = "vox/side" if fam == "absolute" else "% of hemi CC"
        print(f"\n{comp}  {fam} ({unit})   pHGG n={s.n_ped.iloc[0]}")
        print(f"  {'thr':>7} {'pHGG':>7} {'GBM':>7} {'OR':>7} {'p':>8} {'kappa':>7}")
        for _, r in s.iterrows():
            mark = " <-- current" if (fam == "absolute" and
                    ((comp in ("WT","TC") and r.threshold == 20) or (comp == "ET" and r.threshold == 5))) else ""
            print(f"  {r.threshold:>7g} {r.ped_pct:>6.1f}% {r.gbm_pct:>6.1f}% {r.odds_ratio:>7.2f} "
                  f"{ps(r.p):>8} {r.kappa:>7.3f}{mark}")
print("\n" + "="*76)
print("SELECTION on reproducibility (kappa), not on p-value")
print("="*76)
for comp in ["WT", "TC", "ET"]:
    s = S[S.compartment == comp]
    best = s.loc[s.kappa.idxmax()]
    sig = s[s.p < .05]
    print(f"\n{comp}: highest kappa = {best.kappa:.3f} at {best.family} {best.threshold:g}"
          f"  (OR {best.odds_ratio:.2f}, p={ps(best.p)})")
    print(f"    thresholds with kappa >= 0.90: {int((s.kappa>=0.90).sum())}/{len(s)}"
          f"   with p<.05: {len(sig)}/{len(s)}")
    if len(sig):
        print(f"    among p<.05, kappa range {sig.kappa.min():.3f} to {sig.kappa.max():.3f}")
print(f"\nWrote {OUTD/'TableS8_crossing_definition_sweep.csv'}")
