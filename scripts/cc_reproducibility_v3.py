"""
cc_reproducibility_v3.py
Registration reproducibility across independent ANTs random seeds, measured on
the JHU atlas with corrected labels.

Reports
  1. Dice between the reference (seed 42) warped segmentation and each repeat,
     per compartment, so registration stability is quantified rather than assumed
  2. Dice of the callosal intersection itself, which is what the study measures
  3. Agreement of the derived binary outcomes (involvement, midline crossing)
     across seeds, including Cohen kappa and per-seed prevalence
  4. Consensus estimates using majority vote across seeds

Run after cc_multiseed_v3.py completes.
"""
import glob, os, sys
from pathlib import Path
import ants, nibabel as nib, numpy as np, pandas as pd
import config as cfg
OUTD = cfg.FIGURES_SUP; OUTD.mkdir(parents=True, exist_ok=True)
A = cfg.ATLAS
MASK = {r: ants.image_read(str(A/f"CC_{'mask' if r=='whole' else r}.nii.gz")).numpy() > 0
        for r in ["whole", "genu", "body", "splenium"]}
MID = MASK["whole"].shape[0]//2
COH = open(cfg.COHORT).read().split()
SEEDS = [43, 44, 45]
THR = {"WT": 20, "TC": 20, "ET": 5}


def comps(seg):
    return {"WT": seg > 0, "TC": np.isin(seg, [1, 2, 3]), "ET": np.isin(seg, [1])}


def dice(a, b):
    s = a.sum() + b.sum()
    return float(2*np.logical_and(a, b).sum()/s) if s else np.nan


def load(sid, seed):
    p = (cfg.PEDS_REG/f"{sid}_seg_MNI.nii.gz") if seed == 42 \
        else (cfg.seed_reg(seed)/f"{sid}_seg_MNI.nii.gz")
    return np.asanyarray(nib.load(p).dataobj) if p.exists() else None


rows, flags = [], []
for sid in COH:
    ref = load(sid, 42)
    if ref is None:
        continue
    rc = comps(ref)
    rec = {"subject": sid}
    for c, arr in rc.items():
        cc = arr & MASK["whole"]
        rec[f"seed42_{c}_inv"] = int(cc.sum() > 0)
        rec[f"seed42_{c}_cross"] = int(cc[:MID].sum() > THR[c] and cc[MID:].sum() > THR[c])
    for s in SEEDS:
        m = load(sid, s)
        if m is None:
            continue
        mc = comps(m)
        for c in ["WT", "TC", "ET"]:
            rec[f"dice_{c}_s{s}"] = dice(rc[c], mc[c])
            rec[f"diceCC_{c}_s{s}"] = dice(rc[c] & MASK["whole"], mc[c] & MASK["whole"])
            cc = mc[c] & MASK["whole"]
            rec[f"s{s}_{c}_inv"] = int(cc.sum() > 0)
            rec[f"s{s}_{c}_cross"] = int(cc[:MID].sum() > THR[c] and cc[MID:].sum() > THR[c])
    rows.append(rec)

D = pd.DataFrame(rows)
D.to_csv(OUTD/"TableS4_reproducibility_per_case.csv", index=False)
n_full = D[[f"dice_WT_s{s}" for s in SEEDS]].notna().all(axis=1).sum()
print(f"cases with reference + all {len(SEEDS)} repeat registrations: {n_full}/{len(D)}")

print("\n" + "="*86)
print("REGISTRATION REPRODUCIBILITY: Dice vs reference registration (seed 42)")
print("="*86)
print(f"{'compartment':12} {'scope':16} {'median':>8} {'IQR':>18} {'min':>8} {'n':>6}")
summ = []
for c in ["WT", "TC", "ET"]:
    for scope, pref in [("whole brain", "dice"), ("within CC mask", "diceCC")]:
        v = pd.concat([D[f"{pref}_{c}_s{s}"] for s in SEEDS]).dropna()
        if not len(v):
            continue
        print(f"{c:12} {scope:16} {v.median():>8.4f} "
              f"{'['+format(v.quantile(.25),'.4f')+', '+format(v.quantile(.75),'.4f')+']':>18} "
              f"{v.min():>8.4f} {len(v):>6}")
        summ.append(dict(compartment=c, scope=scope, median=v.median(),
                         q1=v.quantile(.25), q3=v.quantile(.75), minimum=v.min(), n=len(v)))
pd.DataFrame(summ).to_csv(OUTD/"TableS5_reproducibility_summary.csv", index=False)

print("\n" + "="*86)
print("OUTCOME STABILITY across seeds (prevalence per seed, and agreement with reference)")
print("="*86)
print(f"{'outcome':22} {'s42':>7} {'s43':>7} {'s44':>7} {'s45':>7} {'consensus':>10} {'agree':>8} {'kappa':>7}")
stab = []
for c in ["WT", "TC", "ET"]:
    for kind in ["inv", "cross"]:
        cols = [f"s{s}_{c}_{kind}" for s in SEEDS if f"s{s}_{c}_{kind}" in D]
        if not cols:
            continue
        sub = D[[f"seed42_{c}_{kind}"]+cols].dropna()
        if not len(sub):
            continue
        ref = sub[f"seed42_{c}_{kind}"].astype(int)
        prev = [100*ref.mean()] + [100*sub[cl].astype(int).mean() for cl in cols]
        allv = pd.concat([ref]+[sub[cl].astype(int) for cl in cols], axis=1)
        cons = (allv.sum(axis=1) >= (allv.shape[1]/2.0))
        agree = np.mean([(ref == sub[cl].astype(int)).mean() for cl in cols])
        ks = []
        for cl in cols:
            a, b = ref.values, sub[cl].astype(int).values
            po = (a == b).mean()
            pe = a.mean()*b.mean() + (1-a.mean())*(1-b.mean())
            ks.append((po-pe)/(1-pe) if pe < 1 else 1.0)
        nm = f"{c} {'involvement' if kind=='inv' else 'crossing'}"
        print(f"{nm:22} " + " ".join(f"{p:>6.1f}%" for p in prev) +
              f" {100*cons.mean():>9.1f}% {100*agree:>7.1f}% {np.mean(ks):>7.3f}")
        stab.append(dict(outcome=nm, prev_s42=prev[0], prev_s43=prev[1] if len(prev)>1 else np.nan,
                         prev_s44=prev[2] if len(prev)>2 else np.nan,
                         prev_s45=prev[3] if len(prev)>3 else np.nan,
                         consensus=100*cons.mean(), mean_agreement=100*agree, mean_kappa=np.mean(ks), n=len(sub)))
pd.DataFrame(stab).to_csv(OUTD/"TableS6_outcome_stability.csv", index=False)
print(f"\nWrote TableS4/S5/S6 to {OUTD}")
