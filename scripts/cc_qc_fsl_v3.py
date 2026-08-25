"""
cc_qc_fsl_v3.py
FSLeyes-style orthographic QC images, one per case, with quantitative overlap.

Each figure: axial / coronal / sagittal through the callosal centroid, JHU
genu/body/splenium contours, tumor compartment overlays, and a printed panel of
Dice and per-label volumes so a mask or label-convention mismatch shows up as a
number rather than requiring visual inspection.
"""
import sys
from pathlib import Path
import ants, nibabel as nib, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import config as cfg
OUT = cfg.FIGURES_SUP; OUT.mkdir(parents=True, exist_ok=True)
A = cfg.ATLAS
CC = ants.image_read(str(A/"CC_mask.nii.gz")).numpy() > 0
G  = ants.image_read(str(A/"CC_genu.nii.gz")).numpy() > 0
B  = ants.image_read(str(A/"CC_body.nii.gz")).numpy() > 0
S  = ants.image_read(str(A/"CC_splenium.nii.gz")).numpy() > 0
T  = nib.load(ants.get_ants_data("mni")).get_fdata()
MID = CC.shape[0]//2
v3 = pd.read_csv(cfg.SIDES).set_index("subject")

def dice(a, b):
    s = a.sum() + b.sum()
    return 2*np.logical_and(a, b).sum()/s if s else float("nan")

def show(ax, bg, overlays, contours, title):
    ax.imshow(bg.T, cmap="gray", origin="lower", aspect="equal")
    for m, c, al in overlays:
        if m.any():
            rgba = np.zeros((*m.shape, 4)); rgba[..., :3] = matplotlib.colors.to_rgb(c); rgba[..., 3] = m*al
            ax.imshow(np.transpose(rgba, (1, 2, 0)).transpose(1, 0, 2) if False else np.transpose(rgba, (1, 0, 2)),
                      origin="lower", aspect="equal")
    for m, c in contours:
        if m.any(): ax.contour(m.T, levels=[.5], colors=c, linewidths=1.2)
    ax.set_title(title, color="white", fontsize=9.5)
    ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_visible(False)

def qc(sid, cohort, idx):
    d = cfg.PEDS_REG if cohort == "pHGG" else cfg.GLI_REG
    seg = np.asanyarray(nib.load(d/f"{sid}_seg_MNI.nii.gz").dataobj)
    if cohort == "pHGG":
        ET, TC = np.isin(seg, [1]), np.isin(seg, [1, 2, 3])
    else:
        ET, TC = np.isin(seg, [3]), np.isin(seg, [1, 3])
    WT = seg > 0
    r = v3.loc[sid]
    inCC = WT & CC
    cx, cy, cz = (int(np.median(np.where(inCC)[i])) for i in range(3)) if inCC.any() \
                 else (MID, CC.shape[1]//2, CC.shape[2]//2)
    fig = plt.figure(figsize=(15.5, 6.6), facecolor="black")
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, .82], wspace=.04)
    ov = lambda sl: [(WT[sl], "#4d8fe8", .42), (TC[sl], "#f4a261", .55), (ET[sl], "#e03a4e", .85)]
    ct = lambda sl: [(G[sl], "#ff3b30"), (B[sl], "#34c759"), (S[sl], "#0a84ff")]
    show(fig.add_subplot(gs[0]), T[:, :, cz], ov(np.s_[:, :, cz]), ct(np.s_[:, :, cz]), f"axial  z={cz-72} mm")
    show(fig.add_subplot(gs[1]), T[:, cy, :], ov(np.s_[:, cy, :]), ct(np.s_[:, cy, :]), f"coronal  y={cy-126} mm")
    ax = fig.add_subplot(gs[2]); show(ax, T[cx, :, :], ov(np.s_[cx, :, :]), ct(np.s_[cx, :, :]), f"sagittal  x={90-cx} mm")
    ax.set_xlim(45, 190); ax.set_ylim(45, 132)
    tx = fig.add_subplot(gs[3]); tx.axis("off"); tx.set_facecolor("black")
    lab = "labels 1=ET 2=NET 3=CC 4=ED" if cohort == "pHGG" else "labels 1=NCR 2=ED 3=ET"
    txt = (f"{sid}\n{cohort}\n\n"
           f"{lab}\n"
           f"  lab1 {int(r.lab1_vol):>7}\n  lab2 {int(r.lab2_vol):>7}\n"
           f"  lab3 {int(r.lab3_vol):>7}\n  lab4 {int(r.lab4_vol):>7}\n\n"
           f"volumes (voxels)\n  WT {int(r.WT_vol):>9}\n  TC {int(r.TC_vol):>9}\n  ET {int(r.ET_vol):>9}\n\n"
           f"DICE with CC atlas\n"
           f"  WT {dice(WT, CC):.4f}\n  TC {dice(TC, CC):.4f}\n  ET {dice(ET, CC):.4f}\n\n"
           f"CC overlap (voxels)\n"
           f"  whole    WT {int(r.whole_WT_n):>6}\n"
           f"  genu     WT {int(r.genu_WT_n):>6}\n"
           f"  body     WT {int(r.body_WT_n):>6}\n"
           f"  splenium WT {int(r.splenium_WT_n):>6}\n\n"
           f"midline crossing\n"
           f"  TC  L={int(r.whole_TC_L):>5} R={int(r.whole_TC_R):>5}\n"
           f"      >20/side: {'YES' if r.whole_TC_L>20 and r.whole_TC_R>20 else 'no'}\n"
           f"  ET  L={int(r.whole_ET_L):>5} R={int(r.whole_ET_R):>5}\n"
           f"      >5/side:  {'YES' if r.whole_ET_L>5 and r.whole_ET_R>5 else 'no'}")
    tx.text(0, 1, txt, color="white", fontsize=8.6, va="top", ha="left", family="monospace",
            transform=tx.transAxes)
    fig.legend(handles=[Patch(color="#4d8fe8", label="whole tumor"), Patch(color="#f4a261", label="tumor core"),
                        Patch(color="#e03a4e", label="enhancing tumor"),
                        Patch(facecolor="none", edgecolor="#ff3b30", label="genu"),
                        Patch(facecolor="none", edgecolor="#34c759", label="body"),
                        Patch(facecolor="none", edgecolor="#0a84ff", label="splenium")],
               loc="lower center", ncol=6, facecolor="black", labelcolor="white",
               fontsize=9, framealpha=0, bbox_to_anchor=(.44, .005))
    fig.suptitle(f"QC {idx}: {sid} ({cohort})", color="white", fontsize=12.5, y=.97)
    fig.subplots_adjust(top=.90, bottom=.09)
    p = OUT/f"QC{idx}_{sid}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="black"); plt.close(fig)
    print("wrote", p)

if __name__ == "__main__":
    ped = v3[v3.cohort == "pHGG"]
    # deliberately span the range: highest involvement, a crossing case,
    # a non-enhancing case, a CC-spared case, and one adult GBM comparator
    hi = ped.sort_values("whole_WT_n", ascending=False).index[0]
    cross = ped[(ped.whole_TC_L > 20) & (ped.whole_TC_R > 20)].sort_values("whole_TC_n", ascending=False).index[0]
    nonenh = ped[ped.ET_vol == 0].sort_values("whole_WT_n", ascending=False).index[0]
    spared = ped[ped.whole_WT_n == 0].sort_values("WT_vol", ascending=False).index[0]
    gbm = v3[v3.cohort == "GBM"].sort_values("whole_TC_n", ascending=False).index[0]
    for i, (sid, coh) in enumerate([(hi, "pHGG"), (cross, "pHGG"), (nonenh, "pHGG"),
                                    (spared, "pHGG"), (gbm, "GBM")], 1):
        qc(sid, coh, i)
