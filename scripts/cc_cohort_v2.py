"""
cc_cohort_v2.py
Corrected inclusion filter and cohort definition.

The original filter (cc_registration_pipeline.py get_phgg_ids) read only the
FIRST cbtn_histologies.tsv row per participant. Patients carry up to 198
biospecimen rows, and 11 CBTN cases record "High-grade glioma/astrocytoma
(WHO grade III/IV)" on a later row while the first row reads NA or Other.
Those cases met the stated inclusion criterion but were dropped.

Corrected filter:
  Training split only
  AND ( Source == DFCI-BCH-BWH-PEDs-HGG
        OR (Source == CBTN AND any row has pathology_diagnosis containing
            "High-grade") )
  AND NOT pontine (no row mentioning diffuse intrinsic pontine glioma, and no
      brainstem primary site) - preserves the pre-existing anatomic exclusion.

Writes cc_cohort_v2.txt and refreshes cc_sides.csv for any new subjects.
"""
import collections, csv, glob, json, os
from pathlib import Path
import config as cfg
NULL = {"", "NA", "Not Applicable", "Not Available", "None"}

meta = list(csv.DictReader(open(cfg.require(cfg.PEDS_META, "BraTS-PEDs metadata")), delimiter="\t"))
rows = collections.defaultdict(list)
for r in csv.DictReader(open(cfg.CBTN), delimiter="\t"):
    c = r["cohort_participant_id"].strip()
    if c:
        rows[c].append(r)


def is_pontine(rs):
    """Pontine veto, matching the original exclusion: DIPG by pathology label only.
    primary_site is multi-valued (e.g. "Thalamus;Brain Stem- Midbrain/Tectum") and
    is NOT used, since thalamic pHGG with midbrain extension is eligible."""
    return any("pontine" in r.get("pathology_diagnosis", "").lower() for r in rs)


included, reasons = [], {}
for r in meta:
    if r["BraTS2025_cohort"] != "Training":
        continue
    sid, src = r["BraTS-SubjectID"], r["Source"]
    rs = rows.get(r["MappingID"].strip(), [])
    if src == "DFCI-BCH-BWH-PEDs-HGG":
        included.append(sid); reasons[sid] = "DFCI collection (curated pHGG)"
    elif src == "CBTN":
        hgg = any("High-grade" in x.get("pathology_diagnosis", "") for x in rs)
        if hgg and not is_pontine(rs):
            first = rs[0].get("pathology_diagnosis", "").strip() if rs else ""
            tag = "CBTN high-grade (first row)" if "High-grade" in first \
                  else "CBTN high-grade (RECOVERED: later biospecimen row)"
            included.append(sid); reasons[sid] = tag
        elif hgg and is_pontine(rs):
            reasons[sid] = "EXCLUDED pontine despite high-grade label"

included = sorted(included)
avail = [s for s in included
         if (cfg.PEDS_TRAIN/s/f"{s}-seg.nii.gz").exists()]

print(f"Corrected filter admits {len(included)}; imaging available for {len(avail)}")
print(collections.Counter(reasons[s] for s in avail).most_common())

rec = sorted(s for s in avail if "RECOVERED" in reasons[s])
print(f"\nRECOVERED cases ({len(rec)}):")
for s in rec:
    print(f"  {s}")
excl_pont = sorted(s for s, v in reasons.items() if v.startswith("EXCLUDED"))
print(f"\nHigh-grade but excluded as pontine ({len(excl_pont)}): {[s[-8:-4] for s in excl_pont]}")

(cfg.COHORT).write_text("\n".join(avail) + "\n")
with open(cfg.COHORT_REASONS, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["subject", "inclusion_reason"])
    for s in avail:
        w.writerow([s, reasons[s]])

# which cases already have a warped segmentation?
have = {f.name.replace("_seg_MNI.nii.gz", "") for f in cfg.PEDS_REG.glob("*_seg_MNI.nii.gz")}
missing = [s for s in avail if s not in have]
print(f"\nwarped segmentations present for {len(set(avail) & have)}/{len(avail)}")
if missing:
    print(f"  {len(missing)} still to register; run scripts/cc_register.py --cohort peds")
print(f"Wrote cc_cohort_v2.txt ({len(avail)} subjects)")
