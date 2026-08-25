"""
cc_qc_v3.py
Rebuilt QC that can actually catch the failures the old QC missed:
  - sagittal panel, not axial only
  - mask CONTOURS, never bounding boxes
  - printed Dice of tumor vs atlas and per-label volumes per case
  - per-case atlas overlap so a mask/registration mismatch is quantified
"""
import glob, os
from pathlib import Path
import ants, nibabel as nib, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import config as cfg
OUT = cfg.FIGURES_SUP; OUT.mkdir(parents=True, exist_ok=True)
A = cfg.ATLAS
CC = ants.image_read(str(A/"CC_mask.nii.gz")).numpy()>0
G  = ants.image_read(str(A/"CC_genu.nii.gz")).numpy()>0
B  = ants.image_read(str(A/"CC_body.nii.gz")).numpy()>0
S  = ants.image_read(str(A/"CC_splenium.nii.gz")).numpy()>0
T  = nib.load(ants.get_ants_data("mni")).get_fdata()
MID = CC.shape[0]//2
v3 = pd.read_csv(cfg.SIDES).set_index("subject")

def panel(fig, gs, i, j, bg, ov, ttl, contours):
    ax = fig.add_subplot(gs[i, j])
    ax.imshow(bg.T, cmap="gray", origin="lower", aspect="equal")
    for m, c, a in ov:
        if m.any():
            rgba = np.zeros((*m.shape,4)); rgba[...,:3]=matplotlib.colors.to_rgb(c); rgba[...,3]=m*a
            ax.imshow(np.transpose(rgba,(1,0,2)), origin="lower", aspect="equal")
    for m, c in contours:
        if m.any(): ax.contour(m.T, levels=[.5], colors=c, linewidths=1.1)
    ax.set_title(ttl, color="white", fontsize=8.5); ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_visible(False)
    return ax

def sheet(cases, cohort, fname, title):
    Z=[y+72 for y in (2,12,22)]
    fig=plt.figure(figsize=(4.1*4, 3.9*len(cases)), facecolor="black")
    gs=fig.add_gridspec(len(cases), 4, wspace=.03, hspace=.16)
    for row,sid in enumerate(cases):
        d = cfg.PEDS_REG if cohort=="pHGG" else cfg.GLI_REG
        seg=np.asanyarray(nib.load(d/f"{sid}_seg_MNI.nii.gz").dataobj)
        ET = np.isin(seg,[1]) if cohort=="pHGG" else np.isin(seg,[3])
        TC = np.isin(seg,[1,2]) if cohort=="pHGG" else np.isin(seg,[1,3])
        WT = seg>0
        r=v3.loc[sid]
        dice = 2*(WT&CC).sum()/(WT.sum()+CC.sum())
        for col,k in enumerate(Z):
            panel(fig,gs,row,col,T[:,:,k],[(WT[:,:,k],"#4d8fe8",.45),(TC[:,:,k],"#f4a261",.6),(ET[:,:,k],"#e03a4e",.85)],
                  f"axial z={k-72}mm",[(G[:,:,k],"#ff3b30"),(B[:,:,k],"#34c759"),(S[:,:,k],"#0a84ff")])
        ax=panel(fig,gs,row,3,T[90,:,:],[(WT[90],"#4d8fe8",.45),(TC[90],"#f4a261",.6),(ET[90],"#e03a4e",.85)],
                 "SAGITTAL x=90",[(G[90],"#ff3b30"),(B[90],"#34c759"),(S[90],"#0a84ff")])
        ax.set_xlim(45,190); ax.set_ylim(45,132)
        lab=(f"{sid}\nDice(WT,CC)={dice:.3f}\n"
             f"labels {r.labels_present}\n"
             f"vols 1:{int(r.lab1_vol)} 2:{int(r.lab2_vol)} 3:{int(r.lab3_vol)} 4:{int(r.lab4_vol)}\n"
             f"CC: WT={int(r.whole_WT_n)} TC={int(r.whole_TC_n)} ET={int(r.whole_ET_n)}\n"
             f"TC L={int(r.whole_TC_L)} R={int(r.whole_TC_R)} cross>20:{'YES' if r.whole_TC_L>20 and r.whole_TC_R>20 else 'no'}")
        fig.text(0.005, 1-(row+0.5)/len(cases), lab, color="white", fontsize=7.4, va="center", ha="left", family="monospace")
    fig.legend(handles=[Patch(color="#4d8fe8",label="whole tumor"),Patch(color="#f4a261",label="tumor core"),
                        Patch(color="#e03a4e",label="enhancing tumor"),
                        Patch(facecolor="none",edgecolor="#ff3b30",label="genu"),
                        Patch(facecolor="none",edgecolor="#34c759",label="body"),
                        Patch(facecolor="none",edgecolor="#0a84ff",label="splenium")],
               loc="lower center", ncol=6, facecolor="black", labelcolor="white", fontsize=9, framealpha=0,
               bbox_to_anchor=(.5,.005))
    fig.suptitle(title, color="white", fontsize=12.5, y=.995)
    fig.subplots_adjust(left=.115, right=.995, top=.955, bottom=.045)
    fig.savefig(OUT/fname, dpi=150, bbox_inches="tight", facecolor="black"); plt.close(fig)
    print("wrote", OUT/fname)

ped=v3[v3.cohort=="pHGG"].sort_values("whole_WT_n",ascending=False)
gbm=v3[v3.cohort=="GBM"].sort_values("whole_WT_n",ascending=False)
sheet(list(ped.index[:5]),"pHGG","FigureS1_qc_phgg.png",
      "Registration and mask QC: pediatric HGG, 5 highest callosal involvement\ncontours = JHU genu/body/splenium; sagittal panel included")
sheet(list(gbm.index[:5]),"GBM","FigureS2_qc_gbm.png",
      "Registration and mask QC: adult GBM, 5 highest callosal involvement\ncontours = JHU genu/body/splenium; sagittal panel included")

# cohort-level QC table
rows=[]
for sid,r in v3.iterrows():
    rows.append(dict(subject=sid, cohort=r.cohort, labels=r.labels_present,
                     WT_vol=int(r.WT_vol), CC_WT=int(r.whole_WT_n),
                     frac_CC=round(r.whole_WT_n/CC.sum(),4)))
q=pd.DataFrame(rows); q.to_csv(OUT/"TableS1_per_case_qc.csv", index=False)
print("wrote", OUT/"TableS1_per_case_qc.csv")
print(f"\ncohort-level check: cases with zero CC overlap  pHGG {int((q[q.cohort=='pHGG'].CC_WT==0).sum())}/114, "
      f"GBM {int((q[q.cohort=='GBM'].CC_WT==0).sum())}/1213")
