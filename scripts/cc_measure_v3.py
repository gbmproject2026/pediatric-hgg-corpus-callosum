"""
cc_measure_v3.py
Measure tumor overlap with the corpus callosum for both cohorts.

Atlas: atlas_jhu/, from the JHU ICBM-DTI-81 white matter atlas (labels 3 genu,
4 body, 5 splenium). See atlas_jhu/PROVENANCE.md.

Segmentation labels differ between the two datasets:
  BraTS-PEDs : 1=ET (enhancing), 2=NET (non-enhancing), 3=CC (cyst), 4=ED (edema)
  BraTS-GLI  : 1=NCR (necrotic), 2=ED (edema), 3=ET (enhancing)

Compartments (BraTS convention):
     WT whole tumor = every labeled voxel
     TC tumor core  = all labeled tumor except peritumoral edema
                      PEDs {1,2,3}   GLI {1,3}
     ET enhancing   = PEDs {1}       GLI {3}

Validation gates abort the run if the atlas or the label conventions are
inconsistent with the data.

Outputs results/cc_sides_v3.csv: per case, per CC region, per compartment, left
and right hemisphere voxel counts plus totals.
"""
import glob, json, os, sys
from multiprocessing import Pool, cpu_count
from pathlib import Path

import ants
import nibabel as nib
import numpy as np
import pandas as pd
import config as cfg
ATLAS = cfg.ATLAS
OUT = cfg.SIDES

REGIONS = ["whole", "genu", "body", "splenium"]
MASK = {r: (ants.image_read(str(ATLAS / f"CC_{'mask' if r=='whole' else r}.nii.gz")).numpy() > 0)
        for r in REGIONS}
SHAPE = MASK["whole"].shape
MID = SHAPE[0] // 2

# BraTS-PEDs 2023 (Kazerooni et al.): four annotated subregions
#   1 = ET  enhancing tumor      (bright on T1CE vs T1 pre-contrast)
#   2 = NET non-enhancing tumor
#   3 = CC  cystic component     (bright T2, dark T1CE)
#   4 = ED  peritumoral edema
# There is NO separate necrotic-core label in BraTS-PEDs.
# The challenge defines pediatric TC as enhancing tumor / cystic component /
# necrosis, i.e. all non-edema tumor. For comparability with the adult TC
# (necrotic core + enhancing tumor, also all non-edema tumor), TC is taken as
# every labeled voxel except peritumoral edema in both cohorts.
LABELS = {
    "pHGG": {"WT": None, "TC": (1, 2, 3), "ET": (1,)},   # WT None = seg>0
    "GBM":  {"WT": None, "TC": (1, 3),    "ET": (3,)},
}

# ---- validation gate 1: atlas sanity -----------------------------------------
def validate_atlas():
    T = nib.load(ants.get_ants_data("mni")).get_fdata()
    brain = T > T.max() * 0.12
    problems = []
    for r in REGIONS:
        m = MASK[r]
        assert m.shape == T.shape, f"{r}: grid mismatch"
        idx = np.array(np.where(m))
        fill = m.sum() / np.prod(idx.max(axis=1) - idx.min(axis=1) + 1)
        ratio = T[m].mean() / T[brain].mean()
        if fill > 0.9:
            problems.append(f"{r}: fill {fill:.2f}, looks like a box")
        if ratio < 1.0:
            problems.append(f"{r}: T1 ratio {ratio:.2f}, not white matter")
    # anterior-to-posterior ordering must be genu > body > splenium
    cy = {r: np.array(np.where(MASK[r]))[1].mean() for r in ["genu", "body", "splenium"]}
    if not (cy["genu"] > cy["body"] > cy["splenium"]):
        problems.append(f"subregion AP ordering wrong: {cy}")
    # subregions must partition the whole mask
    union = MASK["genu"] | MASK["body"] | MASK["splenium"]
    if not np.array_equal(union, MASK["whole"]):
        problems.append("genu|body|splenium != whole mask")
    if MASK["genu"].sum() + MASK["body"].sum() + MASK["splenium"].sum() != MASK["whole"].sum():
        problems.append("subregions overlap each other")
    return problems, cy


def measure(job):
    sid, cohort, path = job
    seg = np.asanyarray(nib.load(path).dataobj)
    if seg.shape != SHAPE:
        return {"subject": sid, "cohort": cohort, "error": f"shape {seg.shape}"}
    row = {"subject": sid, "cohort": cohort, "error": None}
    present = sorted(int(v) for v in np.unique(seg) if v > 0)
    row["labels_present"] = "|".join(map(str, present))
    for l in range(1, 5):
        row[f"lab{l}_vol"] = int((seg == l).sum())
    comp = {}
    for cname, labs in LABELS[cohort].items():
        comp[cname] = (seg > 0) if labs is None else np.isin(seg, labs)
        row[f"{cname}_vol"] = int(comp[cname].sum())
    for r in REGIONS:
        m = MASK[r]
        for cname, arr in comp.items():
            c = arr & m
            row[f"{r}_{cname}_L"] = int(c[:MID].sum())
            row[f"{r}_{cname}_R"] = int(c[MID:].sum())
            row[f"{r}_{cname}_n"] = int(c.sum())
    return row


def build_jobs():
    coh = set(open(cfg.COHORT).read().split())
    jobs = []
    for sid in sorted(coh):
        p = cfg.PEDS_REG / f"{sid}_seg_MNI.nii.gz"
        if p.exists():
            jobs.append((sid, "pHGG", str(p)))
    for f in sorted(glob.glob(str(cfg.GLI_REG/"*_seg_MNI.nii.gz"))):
        jobs.append((os.path.basename(f).replace("_seg_MNI.nii.gz", ""), "GBM", f))
    return jobs, coh


if __name__ == "__main__":
    print("=" * 92)
    print("VALIDATION GATE 1: atlas")
    print("=" * 92)
    probs, cy = validate_atlas()
    for r in REGIONS:
        m = MASK[r]; idx = np.array(np.where(m))
        fill = m.sum() / np.prod(idx.max(axis=1) - idx.min(axis=1) + 1)
        print(f"  {r:9} vox={int(m.sum()):6}  fill={fill:.2f}  centroid y-idx={idx[1].mean():6.1f}")
    print(f"  AP ordering genu({cy['genu']:.0f}) > body({cy['body']:.0f}) > splenium({cy['splenium']:.0f})")
    print(f"  subregions partition whole mask exactly: "
          f"{np.array_equal(MASK['genu']|MASK['body']|MASK['splenium'], MASK['whole'])}")
    if probs:
        print("\n  FAILED:"); [print("   -", p) for p in probs]; sys.exit(1)
    print("  PASS\n")

    jobs, coh = build_jobs()
    npd = sum(1 for j in jobs if j[1] == "pHGG")
    print(f"jobs: {len(jobs)} (pHGG {npd}, GBM {len(jobs)-npd})")
    assert npd == len(coh) == 114, f"cohort coverage {npd} vs {len(coh)}"

    with Pool(max(1, cpu_count() - 2)) as pool:
        rows = pool.map(measure, jobs, chunksize=8)

    df = pd.DataFrame(rows)
    print("\n" + "=" * 92)
    print("VALIDATION GATE 2: measurement integrity")
    print("=" * 92)
    assert df.error.isna().all(), df[df.error.notna()][["subject", "error"]]
    print(f"  no errors across {len(df)} cases")

    ped, gbm = df[df.cohort == "pHGG"], df[df.cohort == "GBM"]
    # label conventions must differ as established
    assert (gbm.lab4_vol == 0).all(), "GBM should have no label 4"
    print(f"  GBM cases with label 4 (should be 0): {int((gbm.lab4_vol>0).sum())}")
    print(f"  pHGG cases with label 4 (edema):      {int((ped.lab4_vol>0).sum())}/{len(ped)}")
    # compartment nesting must hold
    for nm, d in [("pHGG", ped), ("GBM", gbm)]:
        assert (d.ET_vol <= d.TC_vol).all(), f"{nm}: ET > TC"
        assert (d.TC_vol <= d.WT_vol).all(), f"{nm}: TC > WT"
        for r in REGIONS:
            assert (d[f"{r}_ET_n"] <= d[f"{r}_TC_n"]).all(), f"{nm} {r}: ET>TC in region"
            assert (d[f"{r}_TC_n"] <= d[f"{r}_WT_n"]).all(), f"{nm} {r}: TC>WT in region"
        # subregion counts must sum to whole
        for c in ["WT", "TC", "ET"]:
            s = d[f"genu_{c}_n"] + d[f"body_{c}_n"] + d[f"splenium_{c}_n"]
            assert (s == d[f"whole_{c}_n"]).all(), f"{nm} {c}: subregions do not sum to whole"
        # L+R must equal total
        for r in REGIONS:
            for c in ["WT", "TC", "ET"]:
                assert (d[f"{r}_{c}_L"] + d[f"{r}_{c}_R"] == d[f"{r}_{c}_n"]).all(), f"{nm} {r} {c}: L+R != n"
    print("  compartment nesting ET <= TC <= WT: PASS")
    print("  subregions sum to whole mask:       PASS")
    print("  left + right = total:              PASS")

    # cross-check against raw unwarped segmentations for a sample
    T = str(cfg.PEDS_TRAIN)
    chk = 0
    for sid in ped.subject.head(12):
        f = f"{T}/{sid}/{sid}-seg.nii.gz"
        if not os.path.exists(f):
            continue
        raw = np.asanyarray(nib.load(f).dataobj)
        got = set(int(v) for v in np.unique(raw) if v > 0)
        want = set(int(x) for x in ped[ped.subject == sid].labels_present.iloc[0].split("|"))
        assert got == want, f"{sid}: raw {got} vs warped {want}"
        chk += 1
    print(f"  label sets match raw pre-registration segs in {chk}/{chk} sampled cases: PASS")

    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}  ({len(df)} rows, {len(df.columns)} columns)")
    print("\nCohort composition with CORRECTED labels:")
    print(f"  pHGG n={len(ped)}   enhancing (ET>0): {int((ped.ET_vol>0).sum())} "
          f"({100*(ped.ET_vol>0).mean():.1f}%)   non-enhancing: {int((ped.ET_vol==0).sum())}")
    print(f"  GBM  n={len(gbm)}   enhancing (ET>0): {int((gbm.ET_vol>0).sum())} "
          f"({100*(gbm.ET_vol>0).mean():.1f}%)")
    print(f"  median WT volume: pHGG {ped.WT_vol.median():.0f}, GBM {gbm.WT_vol.median():.0f}")
