"""
cc_multiseed_v3.py
Registration reproducibility: re-register the 114-case pediatric cohort with
additional ANTs random seeds and save the warped segmentations so overlap can be
recomputed on the JHU atlas with the corrected labels.

Usage: python cc_multiseed_v3.py 43 44 45
Writes cc_results_seed<N>_v3/registrations/*_seg_MNI.nii.gz
"""
import os, sys, time, json
from multiprocessing import Pool, cpu_count
import ants, numpy as np
import config as cfg
TRAIN = str(cfg.PEDS_TRAIN)
MNI = ants.image_read(ants.get_ants_data("mni"))
COH = open(str(cfg.COHORT)).read().split()


def run(args):
    sid, seed = args
    outdir = str(cfg.seed_reg(seed))
    os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/{sid}_seg_MNI.nii.gz"
    if os.path.exists(out):
        return (sid, seed, "cached")
    t1c, seg = f"{TRAIN}/{sid}/{sid}-t1c.nii.gz", f"{TRAIN}/{sid}/{sid}-seg.nii.gz"
    if not (os.path.exists(t1c) and os.path.exists(seg)):
        return (sid, seed, "missing")
    try:
        t0 = time.time()
        ti, si = ants.image_read(t1c), ants.image_read(seg)
        tb = ants.threshold_image(si, low_thresh=0.5, high_thresh=4.5, inval=1, outval=0)
        td = ants.morphology(tb, operation="dilate", radius=5, mtype="binary", shape="ball")
        cm = ants.threshold_image(td, low_thresh=0.5, inval=0, outval=1)
        reg = ants.registration(fixed=MNI, moving=ti, type_of_transform="SyNRA", mask=cm,
                                reg_iterations=(40, 20, 10), random_seed=seed, verbose=False)
        sm = ants.apply_transforms(fixed=MNI, moving=si, transformlist=reg["fwdtransforms"],
                                   interpolator="nearestNeighbor")
        ants.image_write(sm, out)
        return (sid, seed, f"{time.time()-t0:.1f}s")
    except Exception as e:
        return (sid, seed, f"ERROR {e}")


if __name__ == "__main__":
    seeds = [int(s) for s in sys.argv[1:]] or [43, 44, 45]
    jobs = [(sid, s) for s in seeds for sid in COH]
    print(f"{len(jobs)} registrations ({len(COH)} cases x {len(seeds)} seeds)", flush=True)
    with Pool(max(1, cpu_count() - 2)) as pool:
        res = pool.map(run, jobs, chunksize=4)
    bad = [r for r in res if r[2].startswith(("ERROR", "missing"))]
    print(f"done. errors/missing: {len(bad)}")
    for b in bad[:10]:
        print("  ", b)
