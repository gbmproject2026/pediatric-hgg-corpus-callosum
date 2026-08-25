"""
cc_qc_samples.py
Registration QC panels for the README: sagittal, coronal, axial at the callosal
midline, with the JHU corpus callosum boundary overlaid on the warped tumor.

Cases are chosen to span the range so the computed overlap flags can be checked
against anatomy: CC-involved with midline crossing, involved without crossing,
and CC-spared. Generated for both cohorts.
"""
from pathlib import Path
import ants, nibabel as nib, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import config as cfg
OUT = cfg.QC_SAMPLES
A = cfg.ATLAS
G = ants.image_read(str(A/"CC_genu.nii.gz")).numpy() > 0
B = ants.image_read(str(A/"CC_body.nii.gz")).numpy() > 0
S = ants.image_read(str(A/"CC_splenium.nii.gz")).numpy() > 0
CC = G | B | S
T = nib.load(ants.get_ants_data("mni")).get_fdata()
MID = CC.shape[0]//2
v3 = pd.read_csv(cfg.SIDES).set_index("subject")

def panel(sid, cohort):
    d = cfg.PEDS_REG if cohort == "pHGG" else cfg.GLI_REG
    seg = np.asanyarray(nib.load(d/f"{sid}_seg_MNI.nii.gz").dataobj)
    WT = seg > 0
    TC = np.isin(seg, [1, 2, 3]) if cohort == "pHGG" else np.isin(seg, [1, 3])
    r = v3.loc[sid]
    inCC = WT & CC
    if inCC.any():
        cx, cy, cz = (int(np.median(np.where(inCC)[i])) for i in range(3))
    else:
        cx, cy, cz = MID, 120, 90
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), facecolor="black")
    views = [("sagittal", T[cx, :, :], np.s_[cx, :, :], (45, 190), (45, 132)),
             ("coronal",  T[:, cy, :], np.s_[:, cy, :], (35, 148), (35, 150)),
             ("axial",    T[:, :, cz], np.s_[:, :, cz], (35, 148), (30, 200))]
    for ax, (nm, bg, sl, xl, yl) in zip(axes, views):
        ax.imshow(bg.T, cmap="gray", origin="lower", aspect="equal")
        m = WT[sl]
        if m.any():
            rgba = np.zeros((*m.shape, 4)); rgba[..., :3] = (0.86, 0.15, 0.20); rgba[..., 3] = m*0.55
            ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", aspect="equal")
        for M, c in [(G[sl], "#ff375f"), (B[sl], "#30d158"), (S[sl], "#0a84ff")]:
            if M.any(): ax.contour(M.T, levels=[.5], colors=c, linewidths=1.5)
        ax.set_title(nm, color="white", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(*xl); ax.set_ylim(*yl)
        for s_ in ax.spines.values(): s_.set_visible(False)
    cross = "yes" if (r.whole_TC_L > 20 and r.whole_TC_R > 20) else "no"
    fig.suptitle(f"{sid}  ({cohort})    CC whole (WT) {int(r.whole_WT_n)} vox    "
                 f"genu {int(r.genu_WT_n)}  body {int(r.body_WT_n)}  splenium {int(r.splenium_WT_n)}    "
                 f"TC midline crossing: {cross}", color="white", fontsize=11.5, y=0.97)
    fig.legend(handles=[Patch(color="#dc2626", label="warped tumor (WT)"),
                        Patch(facecolor="none", edgecolor="#ff375f", label="genu"),
                        Patch(facecolor="none", edgecolor="#30d158", label="body"),
                        Patch(facecolor="none", edgecolor="#0a84ff", label="splenium")],
               loc="lower center", ncol=4, facecolor="black", labelcolor="white",
               fontsize=9.5, framealpha=0, bbox_to_anchor=(.5, .01))
    fig.subplots_adjust(top=.88, bottom=.10, wspace=.02)
    p = OUT/f"QC_{sid}.png"
    fig.savefig(p, dpi=130, bbox_inches="tight", facecolor="black"); plt.close(fig)
    return p, int(r.whole_WT_n), int(r.genu_WT_n), int(r.splenium_WT_n), cross

if __name__ == "__main__":
    ped = v3[v3.cohort == "pHGG"]; gbm = v3[v3.cohort == "GBM"]
    pe = ped[ped.ET_vol > 0]
    picks = []
    picks.append((pe[(pe.whole_TC_L > 20) & (pe.whole_TC_R > 20)].sort_values("splenium_WT_n", ascending=False).index[0], "pHGG"))
    picks.append((ped[(ped.genu_WT_n > 0)].sort_values("genu_WT_n", ascending=False).index[0], "pHGG"))
    picks.append((ped[(ped.whole_WT_n > 0) & (ped.whole_TC_L <= 20)].sort_values("whole_WT_n", ascending=False).index[0], "pHGG"))
    picks.append((ped[ped.whole_WT_n == 0].sort_values("WT_vol", ascending=False).index[0], "pHGG"))
    picks.append((gbm[(gbm.whole_TC_L > 20) & (gbm.whole_TC_R > 20)].sort_values("genu_WT_n", ascending=False).index[0], "GBM"))
    picks.append((gbm[gbm.whole_WT_n == 0].sort_values("WT_vol", ascending=False).index[0], "GBM"))
    print(f"{'case':24} {'cohort':7} {'CC WT':>7} {'genu':>7} {'splen':>7} {'cross':>6}")
    for sid, coh in picks:
        p, w, g, s, c = panel(sid, coh)
        print(f"{sid:24} {coh:7} {w:>7} {g:>7} {s:>7} {c:>6}")
    print(f"\nwrote {len(picks)} panels to {OUT}")
