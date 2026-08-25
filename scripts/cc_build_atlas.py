"""
cc_build_atlas.py
Build atlas_jhu/ from the JHU ICBM-DTI-81 white matter atlas.

The atlas ships 48 white matter labels; the corpus callosum occupies three of
them (3 genu, 4 body, 5 splenium). It is already on the 182x218x182 MNI152 grid
used by ants.get_ants_data('mni'), so no resampling or interpolation is needed.

Every check below must pass before the masks are written. See
docs/METHODS_NOTES.md for what each one verifies.

Usage: python cc_build_atlas.py [--source PATH_OR_URL] [--out DIR]
"""
import argparse, sys, urllib.request
from pathlib import Path
import ants, nibabel as nib, numpy as np

DEFAULT_SRC = ("https://raw.githubusercontent.com/neurodata/neuroparc/master/"
               "atlases/label/Human/JHU_space-MNI152NLin6_res-1x1x1.nii.gz")
LABELS = {"genu": 3, "body": 4, "splenium": 5}


def fetch(src, dst):
    if str(src).startswith("http"):
        print(f"downloading {src}")
        urllib.request.urlretrieve(src, dst)
    else:
        dst.write_bytes(Path(src).read_bytes())
    return dst


def build(src, out):
    out.mkdir(parents=True, exist_ok=True)
    jhu = out/"JHU_ICBM-DTI-81_space-MNI152NLin6_res-1mm.nii.gz"
    if not jhu.exists():
        fetch(src, jhu)
    img = nib.load(jhu)
    J = np.asanyarray(img.dataobj).astype(int)

    tmpl = nib.load(ants.get_ants_data("mni"))
    assert J.shape == tmpl.shape, f"grid mismatch: atlas {J.shape} vs template {tmpl.shape}"
    assert np.allclose(img.affine, tmpl.affine), "affine mismatch with the registration target"
    print(f"grid matches registration target: {J.shape}")

    T = tmpl.get_fdata(); brain = T > T.max()*0.12
    masks = {}
    for name, lab in LABELS.items():
        m = (J == lab)
        assert m.any(), f"label {lab} ({name}) absent from the atlas"
        masks[name] = m
        nib.save(nib.Nifti1Image(m.astype(np.uint8), img.affine, dtype=np.uint8), out/f"CC_{name}.nii.gz")
    whole = masks["genu"] | masks["body"] | masks["splenium"]
    nib.save(nib.Nifti1Image(whole.astype(np.uint8), img.affine, dtype=np.uint8), out/"CC_mask.nii.gz")

    print(f"\n{'region':10} {'voxels':>8} {'fill':>6} {'T1 ratio':>9} {'centroid y idx':>15}")
    for name, m in list(masks.items()) + [("whole", whole)]:
        idx = np.array(np.where(m))
        fill = m.sum()/np.prod(idx.max(axis=1) - idx.min(axis=1) + 1)
        ratio = T[m].mean()/T[brain].mean()
        print(f"{name:10} {int(m.sum()):>8} {fill:>6.2f} {ratio:>9.2f} {idx[1].mean():>15.1f}")
        # a shaped structure, not a bounding box
        assert fill < 0.9, f"{name}: fill {fill:.2f}, expected a shaped structure"
        # corpus callosum is white matter and must be brighter than whole brain
        assert ratio > 1.0, f"{name}: T1 ratio {ratio:.2f}, not white matter"

    cy = {n: np.array(np.where(masks[n]))[1].mean() for n in LABELS}
    assert cy["genu"] > cy["body"] > cy["splenium"], (
        f"anterior-posterior ordering wrong: {cy}. Higher index is anterior; "
        "the genu is the anterior bend and the splenium the posterior bulb.")
    print(f"\nanterior-to-posterior ordering genu > body > splenium: OK")

    assert masks["genu"].sum() + masks["body"].sum() + masks["splenium"].sum() == whole.sum(), \
        "subregions overlap each other"
    assert np.array_equal(masks["genu"] | masks["body"] | masks["splenium"], whole), \
        "subregions do not partition the whole mask"
    print("subregions partition the whole mask exactly: OK")
    print(f"\nwrote {out}/CC_mask.nii.gz and the three subregions")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=DEFAULT_SRC)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent/"atlas_jhu"))
    a = ap.parse_args()
    build(a.source, Path(a.out))
