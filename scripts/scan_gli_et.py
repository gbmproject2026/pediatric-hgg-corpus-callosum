"""
scan_gli_et.py
Define the adult comparator cohort.

The BraTS-GLI download contains cases of mixed grade. This selects cases with
more than 100 enhancing-tumor voxels, using enhancement as a proxy for
glioblastoma, and writes the surviving subject ids for cc_register.py.

Adult BraTS-GLI labels: 1 necrotic core, 2 peritumoral edema, 3 enhancing tumor.
Enhancing tumor is label 3 here. Do not apply this to BraTS-PEDs, where label 3
is the cystic component; see docs/METHODS_NOTES.md.

Usage:
  python scan_gli_et.py --train-dir /path/BraTS-GLI --out ../results/gbm_cases.txt
"""
import argparse, glob, os
from multiprocessing import Pool, cpu_count
import nibabel as nib, numpy as np

ET_LABEL = 3
MIN_ET_VOXELS = 100


def count(seg_path):
    sid = os.path.basename(os.path.dirname(seg_path))
    try:
        a = np.asanyarray(nib.load(seg_path).dataobj)
        return sid, int((a == ET_LABEL).sum())
    except Exception as e:
        return sid, -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-et", type=int, default=MIN_ET_VOXELS)
    a = ap.parse_args()

    segs = sorted(glob.glob(os.path.join(a.train_dir, "*", "*-seg.nii.gz")))
    if not segs:
        segs = sorted(glob.glob(os.path.join(a.train_dir, "*", "*", "*-seg.nii.gz")))
    print(f"scanning {len(segs)} segmentations for enhancing tumor (label {ET_LABEL})")

    with Pool(max(1, cpu_count()-2)) as pool:
        res = pool.map(count, segs, chunksize=8)

    ok = [(s, n) for s, n in res if n > a.min_et]
    zero = sum(1 for _, n in res if n == 0)
    low = sum(1 for _, n in res if 0 < n <= a.min_et)
    bad = sum(1 for _, n in res if n < 0)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        f.write("\n".join(s for s, _ in sorted(ok)) + "\n")

    print(f"\nkept    {len(ok)}  (ET > {a.min_et} voxels)")
    print(f"dropped {zero} with no enhancement, {low} at or below threshold, {bad} unreadable")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
