"""
cc_register.py
Register a cohort to MNI152 and save the warped segmentations.

This script does registration ONLY. Corpus callosum overlap is computed
afterwards by cc_measure_v3.py, which applies the atlas and the per-dataset label
conventions. Keeping the two stages separate means any change to the atlas or the
label mapping is a re-run of measurement alone, not of registration.

Per case:
  1. read the post-contrast T1 and the expert segmentation
  2. build a cost-function mask by dilating the tumor 5 mm, so tumor tissue does
     not bias alignment
  3. register T1c to the ANTs MNI152 template with SyNRA and a fixed random seed
  4. apply the transform to the segmentation with nearest-neighbor interpolation
     so label values are preserved exactly
  5. write <out>/registrations/<subject>_seg_MNI.nii.gz

Existing outputs are skipped, so an interrupted run resumes on re-run.

Usage
  python cc_register.py --cohort peds --train-dir /path/BraTS-PEDs-v1/Training \
      --ids ../results/cc_cohort_v2.txt --out ../cc_results_peds
  python cc_register.py --cohort gbm --train-dir /path/BraTS-GLI \
      --ids ../results/gbm_cases.txt --out ../cc_results_gbm
"""
import argparse, os, sys, time
from multiprocessing import Pool, cpu_count
from pathlib import Path

import ants


def register(job):
    sid, train_dir, outdir, seed = job
    out = Path(outdir)/"registrations"/f"{sid}_seg_MNI.nii.gz"
    if out.exists():
        return (sid, "cached")
    t1c = Path(train_dir)/sid/f"{sid}-t1c.nii.gz"
    seg = Path(train_dir)/sid/f"{sid}-seg.nii.gz"
    if not (t1c.exists() and seg.exists()):
        return (sid, "missing input")
    try:
        t0 = time.time()
        mni = ants.image_read(ants.get_ants_data("mni"))
        ti, si = ants.image_read(str(t1c)), ants.image_read(str(seg))
        tb = ants.threshold_image(si, low_thresh=0.5, high_thresh=4.5, inval=1, outval=0)
        td = ants.morphology(tb, operation="dilate", radius=5, mtype="binary", shape="ball")
        cm = ants.threshold_image(td, low_thresh=0.5, inval=0, outval=1)
        reg = ants.registration(fixed=mni, moving=ti, type_of_transform="SyNRA",
                                mask=cm, reg_iterations=(40, 20, 10),
                                random_seed=seed, verbose=False)
        warped = ants.apply_transforms(fixed=mni, moving=si,
                                       transformlist=reg["fwdtransforms"],
                                       interpolator="nearestNeighbor")
        out.parent.mkdir(parents=True, exist_ok=True)
        ants.image_write(warped, str(out))
        return (sid, f"{time.time()-t0:.1f}s")
    except Exception as e:
        return (sid, f"ERROR {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", choices=["peds", "gbm"], required=True)
    ap.add_argument("--train-dir", required=True,
                    help="directory holding one folder per subject")
    ap.add_argument("--ids", required=True,
                    help="text file with one subject id per line")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=max(1, cpu_count()-2))
    a = ap.parse_args()

    ids = [x.strip() for x in open(a.ids) if x.strip()]
    if not ids:
        sys.exit(f"no subject ids in {a.ids}")
    print(f"{a.cohort}: {len(ids)} cases, seed {a.seed}, {a.workers} workers")
    print(f"  input  {a.train_dir}")
    print(f"  output {a.out}/registrations")

    jobs = [(s, a.train_dir, a.out, a.seed) for s in ids]
    with Pool(a.workers) as pool:
        res = pool.map(register, jobs, chunksize=2)

    cached = sum(1 for _, r in res if r == "cached")
    bad = [(s, r) for s, r in res if r.startswith(("ERROR", "missing"))]
    print(f"\ndone: {len(res)-len(bad)} available ({cached} already present), {len(bad)} failed")
    for s, r in bad[:20]:
        print(f"  {s}: {r}")
    if bad:
        print("\nRe-run to retry failures; successful cases are skipped.")


if __name__ == "__main__":
    main()
