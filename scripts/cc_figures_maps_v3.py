"""cc_figures_maps_v3.py  Figures 1, 5, 6 on the JHU atlas with corrected labels."""
import glob, os
from pathlib import Path
import ants, nibabel as nib, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import ndimage
import config as cfg
OUT = cfg.FIGURES; OUT.mkdir(parents=True, exist_ok=True)
A = cfg.ATLAS
G = ants.image_read(str(A/"CC_genu.nii.gz")).numpy() > 0
B = ants.image_read(str(A/"CC_body.nii.gz")).numpy() > 0
S = ants.image_read(str(A/"CC_splenium.nii.gz")).numpy() > 0
CC = G | B | S
T = nib.load(ants.get_ants_data("mni")).get_fdata()
MID = CC.shape[0]//2
COH = open(cfg.COHORT).read().split()
v3 = pd.read_csv(cfg.SIDES).set_index("subject")

# ---- FIG 1: atlas ----------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))
for ax, (kind, idx, ttl) in zip(axes, [("sag", 90, "midsagittal x=0 mm"),
                                       ("sag", 96, "parasagittal x=-6 mm"),
                                       ("ax", 22+72, "axial z=22 mm")]):
    bg = T[idx, :, :] if kind == "sag" else T[:, :, idx]
    sl = np.s_[idx, :, :] if kind == "sag" else np.s_[:, :, idx]
    ax.imshow(bg.T, cmap="gray", origin="lower", aspect="equal")
    for m, c in [(G, "#e63946"), (B, "#2a9d8f"), (S, "#4361ee")]:
        if m[sl].any():
            rgba = np.zeros((*m[sl].shape, 4)); rgba[..., :3] = matplotlib.colors.to_rgb(c); rgba[..., 3] = m[sl]*.65
            ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", aspect="equal")
    if kind == "sag":
        ax.set_xlim(45, 190); ax.set_ylim(45, 132)
    else:
        ax.set_xlim(52, 130); ax.set_ylim(50, 180)
    ax.set_title(ttl, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
fig.legend(handles=[Patch(color="#e63946", label="genu (JHU label 3)"),
                    Patch(color="#2a9d8f", label="body (JHU label 4)"),
                    Patch(color="#4361ee", label="splenium (JHU label 5)")],
           loc="lower center", ncol=3, fontsize=10.5, frameon=False)
fig.suptitle("Corpus callosum atlas used for measurement (JHU ICBM-DTI-81, MNI152 space)", fontsize=13)
fig.subplots_adjust(bottom=.13, top=.87, wspace=.04)
fig.savefig(OUT/"Figure1_atlas.png", dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
print("wrote Figure1_atlas.png")

# ---- FIG 5: frequency difference -------------------------------------------
NEIGH = ndimage.binary_dilation(CC, iterations=6)
def freq(ids, d):
    acc = np.zeros(CC.shape, np.float32); n = 0
    for sid in ids:
        f = d/f"{sid}_seg_MNI.nii.gz"
        if not f.exists(): continue
        acc += (np.asanyarray(nib.load(f).dataobj) > 0).astype(np.float32); n += 1
    return acc/max(n, 1), n
gli_ids = [os.path.basename(f).replace("_seg_MNI.nii.gz", "")
           for f in glob.glob(str(cfg.GLI_REG/"*_seg_MNI.nii.gz"))]
pf, npd = freq(COH, cfg.PEDS_REG)
gf, ngl = freq(gli_ids, cfg.GLI_REG)
diff = (pf - gf)*100.0; diff[~NEIGH] = np.nan
vmax = np.nanpercentile(np.abs(diff), 99)
Z = [y+72 for y in (2, 12, 22, 30)]
fig, axes = plt.subplots(1, len(Z)+1, figsize=(4.0*(len(Z)+1), 5.2))
for i, k in enumerate(Z):
    ax = axes[i]
    ax.imshow(T[:, :, k].T, cmap="gray", origin="lower", aspect="equal")
    im = ax.imshow(np.ma.masked_invalid(diff[:, :, k]).T, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   origin="lower", aspect="equal", alpha=.92)
    ax.contour(CC[:, :, k].T, levels=[.5], colors="#111111", linewidths=.9)
    ax.set_title(f"axial z={k-72} mm", fontsize=10.5)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(52, 130); ax.set_ylim(50, 180)
ax = axes[-1]
ax.imshow(T[90, :, :].T, cmap="gray", origin="lower", aspect="equal")
ax.imshow(np.ma.masked_invalid(diff[90, :, :]).T, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
          origin="lower", aspect="equal", alpha=.92)
ax.contour(CC[90, :, :].T, levels=[.5], colors="#111111", linewidths=.9)
ax.set_title("midsagittal x=0 mm", fontsize=10.5)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(45, 190); ax.set_ylim(45, 132)
cb = fig.colorbar(im, ax=axes.tolist(), fraction=.016, pad=.01)
cb.set_label("difference in tumor frequency, percentage points", fontsize=9.5)
fig.suptitle(f"Voxelwise tumor frequency, pediatric HGG (n={npd}) minus adult GBM (n={ngl})\n"
             f"warm = higher in pediatric; black contour = corpus callosum", fontsize=12.5, y=1.0)
fig.savefig(OUT/"Figure5_frequency_difference.png", dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
print("wrote Figure5_frequency_difference.png")

# ---- FIG 6: representative cases -------------------------------------------
ped = v3[v3.cohort == "pHGG"]; gbm = v3[v3.cohort == "GBM"]
p_sel = ped[(ped.whole_TC_L > 20) & (ped.whole_TC_R > 20)].sort_values("whole_TC_n", ascending=False).index[0]
g_sel = gbm[(gbm.whole_TC_L > 20) & (gbm.whole_TC_R > 20)].sort_values("whole_TC_n", ascending=False).index[0]
fig, axes = plt.subplots(2, 4, figsize=(16.5, 8.6), facecolor="black")
for row, (sid, coh, d) in enumerate([(p_sel, "Pediatric HGG", cfg.PEDS_REG),
                                     (g_sel, "Adult GBM", cfg.GLI_REG)]):
    seg = np.asanyarray(nib.load(d/f"{sid}_seg_MNI.nii.gz").dataobj)
    if coh.startswith("Ped"):
        ET, TC = np.isin(seg, [1]), np.isin(seg, [1, 2, 3])
    else:
        ET, TC = np.isin(seg, [3]), np.isin(seg, [1, 3])
    WT = seg > 0
    r = v3.loc[sid]
    for col, k in enumerate(Z):
        ax = axes[row, col]
        ax.imshow(T[:, :, k].T, cmap="gray", origin="lower", aspect="equal")
        for m, c, al in [(WT[:, :, k], "#4d8fe8", .40), (TC[:, :, k], "#f4a261", .55), (ET[:, :, k], "#e03a4e", .85)]:
            if m.any():
                rgba = np.zeros((*m.shape, 4)); rgba[..., :3] = matplotlib.colors.to_rgb(c); rgba[..., 3] = m*al
                ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", aspect="equal")
        for m, c in [(G[:, :, k], "#ff3b30"), (B[:, :, k], "#34c759"), (S[:, :, k], "#0a84ff")]:
            if m.any(): ax.contour(m.T, levels=[.5], colors=c, linewidths=1.1)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(52, 130); ax.set_ylim(50, 180)
        if row == 0: ax.set_title(f"z={k-72} mm", color="white", fontsize=10)
        if col == 0:
            ax.text(-.08, .5, f"{coh}\n{sid}\nTC in CC = {int(r.whole_TC_n)} vox\nL={int(r.whole_TC_L)} R={int(r.whole_TC_R)}",
                    transform=ax.transAxes, color="white", fontsize=9, rotation=90, ha="center", va="center")
fig.legend(handles=[Patch(color="#4d8fe8", label="whole tumor"), Patch(color="#f4a261", label="tumor core"),
                    Patch(color="#e03a4e", label="enhancing tumor"),
                    Patch(facecolor="none", edgecolor="#ff3b30", label="genu"),
                    Patch(facecolor="none", edgecolor="#34c759", label="body"),
                    Patch(facecolor="none", edgecolor="#0a84ff", label="splenium")],
           loc="lower center", ncol=6, facecolor="black", labelcolor="white", fontsize=9.5,
           framealpha=0, bbox_to_anchor=(.5, .02))
fig.suptitle("Representative cases with tumor core crossing the midline", color="white", fontsize=13, y=.97)
fig.subplots_adjust(top=.91, bottom=.08, wspace=.03, hspace=.05)
fig.savefig(OUT/"Figure6_representative.png", dpi=200, bbox_inches="tight", facecolor="black"); plt.close(fig)
print("wrote Figure6_representative.png")
