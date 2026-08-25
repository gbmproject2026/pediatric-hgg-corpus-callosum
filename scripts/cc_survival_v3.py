"""
cc_survival_v3.py
Exploratory survival analysis on the corrected n=114 cohort, linked to CBTN
histology, methylation class, MGMT status, and germline predisposition.

Uses the JHU corpus callosum atlas and the corrected BraTS-PEDs label convention
(ET=1, NET=2, CC=3, ED=4; TC = all tumor except edema).

Outputs
  results/cc_survival_v3.csv                     per-patient table
  figures/supplemental/FigureS3_survival.png     KM panels + hazard forest
  figures/supplemental/FigureS4_molecular.png    molecular subgroup survival
"""
import collections, csv, math
from pathlib import Path
import numpy as np, pandas as pd, scipy.stats as sps
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import config as cfg
OUT = cfg.FIGURES_SUP; OUT.mkdir(parents=True, exist_ok=True)
NULL = {"", "NA", "Not Applicable", "Not Available", "None", "unknown", "Unknown"}
THR = {"WT": 20, "TC": 20, "ET": 5}

v3 = pd.read_csv(cfg.SIDES).set_index("subject")
coh = open(cfg.COHORT).read().split()
meta = {r["BraTS-SubjectID"]: r for r in
        csv.DictReader(open(str(cfg.PEDS_META)), delimiter="\t")}
rows_by_pt = collections.defaultdict(list)
for r in csv.DictReader(open(cfg.CBTN), delimiter="\t"):
    c = r["cohort_participant_id"].strip()
    if c: rows_by_pt[c].append(r)

def pick(rs, col):
    for r in rs:
        v = r.get(col, "").strip()
        if v not in NULL: return v
    return None

K27 = {"DMG_K27"}
PHGG_OTHER = {"pedHGG_RTK1A","pedHGG_RTK1B","pedHGG_RTK2A","pedHGG_RTK2B","pedHGG_MYCN",
              "DHG_G34","IHG","GBM_RTK1","GBM_RTK2","GBM_MES_TYP","A_IDH_HG","DMG_EGFR","pedHGG_B","HGAP"}
def mgroup(sc, score):
    if sc is None: return "no call"
    if score is not None and score < 0.9: return "low confidence"
    if sc in K27: return "DMG H3 K27-altered"
    if sc in PHGG_OTHER: return "other pHGG"
    return "outside pHGG (WHO 2021)"

recs = []
for sid in coh:
    rs = rows_by_pt.get(meta.get(sid, {}).get("MappingID", "").strip(), [])
    if not rs: continue
    od, os_ = pick(rs, "OS_days"), pick(rs, "OS_status")
    if od is None or os_ not in ("DECEASED", "LIVING"): continue
    r = v3.loc[sid]
    sc = pick(rs, "dkfz_v12_methylation_subclass")
    scr = pick(rs, "dkfz_v12_methylation_subclass_score")
    scrf = float(scr) if scr else None
    age = pick(rs, "age_at_diagnosis_days")
    efs = pick(rs, "EFS_days")
    rec = dict(subject=sid, os_days=float(od), os_event=int(os_ == "DECEASED"),
               efs_days=float(efs) if efs else np.nan, efs_event=1 if efs else 0,
               age_yrs=round(float(age)/365.25, 1) if age else np.nan,
               sex=pick(rs, "reported_gender"), cns_region=pick(rs, "CNS_region"),
               resection=pick(rs, "extent_of_tumor_resection"),
               histology=pick(rs, "pathology_diagnosis"),
               free_text=(pick(rs, "pathology_free_text_diagnosis") or "")[:70],
               meth_subclass=sc, meth_score=scrf,
               meth_group=mgroup(sc, scrf),
               mgmt=pick(rs, "dkfz_v12_methylation_mgmt_status"),
               predisposition=pick(rs, "cancer_predispositions"),
               WT_vol=int(r.WT_vol), ET_vol=int(r.ET_vol), enhancing=int(r.ET_vol > 0))
    for c in ["WT", "TC", "ET"]:
        rec[f"inv_{c}"] = int(r[f"whole_{c}_n"] > 0)
        rec[f"cross_{c}"] = int(r[f"whole_{c}_L"] > THR[c] and r[f"whole_{c}_R"] > THR[c])
    for reg in ["genu", "body", "splenium"]:
        rec[f"inv_{reg}_WT"] = int(r[f"{reg}_WT_n"] > 0)
        rec[f"inv_{reg}_TC"] = int(r[f"{reg}_TC_n"] > 0)
    rec["cc_frac"] = round(r.whole_WT_n/35291, 4)
    recs.append(rec)

D = pd.DataFrame(recs)
D["log_wt"] = np.log10(D.WT_vol.clip(lower=1))
D.to_csv(cfg.SURVIVAL, index=False)
D.to_csv(OUT/"DataS4_survival_per_patient.csv", index=False)

print("="*96)
print(f"SURVIVAL SUBSET  n={len(D)} of 114   deaths={int(D.os_event.sum())}  censored={int((1-D.os_event).sum())}")
print(f"  median OS {D.os_days.median():.0f} d   median EFS {D.efs_days.median():.0f} d")
print(f"  median age {D.age_yrs.median():.1f} y   enhancing {int(D.enhancing.sum())}/{len(D)}")
print("="*96)
print("\nHistology (pathology_diagnosis):")
for k, v in D.histology.value_counts().items(): print(f"  {v:3}  {k[:70]}")
print("\nMethylation group:")
for k, v in D.meth_group.value_counts().items(): print(f"  {v:3}  {k}")
print("\nMethylation subclass (confident calls, score>=0.9):")
for k, v in D[D.meth_score >= 0.9].meth_subclass.value_counts().items(): print(f"  {v:3}  {k}")
print("\nMGMT:", dict(D.mgmt.value_counts(dropna=False).items()))
print("Germline predisposition:", dict(D.predisposition.value_counts(dropna=False).items()))
print("Resection:", dict(D.resection.value_counts(dropna=False).items()))

def surv_row(d, key, tcol, ecol, label):
    g1, g0 = d[d[key] == 1], d[d[key] == 0]
    if len(g1) < 3 or len(g0) < 3: return None
    lr = logrank_test(g1[tcol], g0[tcol], g1[ecol], g0[ecol])
    k1, k0 = KaplanMeierFitter(), KaplanMeierFitter()
    k1.fit(g1[tcol], g1[ecol]); k0.fit(g0[tcol], g0[ecol])
    c = CoxPHFitter().fit(d[[key, tcol, ecol]], tcol, ecol)
    hr = float(np.exp(c.params_[key])); lo, hi = (float(np.exp(x)) for x in c.confidence_intervals_.loc[key])
    try:
        c2 = CoxPHFitter().fit(d[[key, "log_wt", tcol, ecol]], tcol, ecol)
        ahr = float(np.exp(c2.params_[key])); ap = float(c2.summary.loc[key, "p"])
    except Exception:
        ahr, ap = np.nan, np.nan
    return dict(label=label, n1=len(g1), n0=len(g0), e1=int(g1[ecol].sum()), e0=int(g0[ecol].sum()),
                med1=k1.median_survival_time_, med0=k0.median_survival_time_,
                logrank_p=lr.p_value, hr=hr, lo=lo, hi=hi,
                p=float(c.summary.loc[key, "p"]), ahr=ahr, ap=ap)

PREDS = [("inv_WT","CC involvement (WT)"),("inv_TC","CC involvement (TC)"),
         ("cross_WT","Midline crossing (WT)"),("cross_TC","Midline crossing (TC)"),
         ("inv_genu_WT","Genu involvement (WT)"),("inv_body_WT","Body involvement (WT)"),
         ("inv_splenium_WT","Splenium involvement (WT)")]
def ps(p): return "<.001" if p < .001 else f"{p:.3f}".lstrip("0")

res = {}
for tcol, ecol, nm in [("os_days","os_event","OVERALL SURVIVAL"),("efs_days","efs_event","EVENT-FREE SURVIVAL")]:
    d = D.dropna(subset=[tcol]); d = d[d[tcol] > 0]
    print(f"\n{'='*96}\n{nm}   n={len(d)}, events={int(d[ecol].sum())}\n{'='*96}")
    print(f"{'predictor':26} {'inv/spared':>12} {'median (d)':>16} {'log-rank':>9} {'HR (95% CI)':>22} {'adj HR':>8} {'adj p':>7}")
    out = []
    for k, lab in PREDS:
        r = surv_row(d, k, tcol, ecol, lab)
        if r is None:
            print(f"{lab:26} group too small"); continue
        out.append(r)
        print(f"{lab:26} {str(r['n1'])+'/'+str(r['n0']):>12} "
              f"{f'{r[chr(109)+chr(101)+chr(100)+chr(49)]:.0f} vs {r[chr(109)+chr(101)+chr(100)+chr(48)]:.0f}':>16} "
              f"{ps(r['logrank_p']):>9} {f'{r[chr(104)+chr(114)]:.2f} ({r[chr(108)+chr(111)]:.2f}-{r[chr(104)+chr(105)]:.2f})':>22} "
              f"{r['ahr']:>8.2f} {ps(r['ap']):>7}")
    res[nm] = pd.DataFrame(out)
    res[nm].to_csv(OUT/f"TableS7_survival_{'os' if 'OVERALL' in nm else 'efs'}.csv", index=False)

# molecular
print(f"\n{'='*96}\nMOLECULAR SUBGROUP\n{'='*96}")
conf = D[(D.meth_score >= 0.9) & (D.meth_group.isin(["DMG H3 K27-altered","other pHGG"]))]
g1, g0 = conf[conf.meth_group == "DMG H3 K27-altered"], conf[conf.meth_group == "other pHGG"]
if len(g1) >= 3 and len(g0) >= 3:
    lr = logrank_test(g1.os_days, g0.os_days, g1.os_event, g0.os_event)
    for nm2, g in [("DMG H3 K27-altered", g1), ("other pHGG", g0)]:
        k = KaplanMeierFitter(); k.fit(g.os_days, g.os_event)
        print(f"  {nm2:22} n={len(g):2} deaths={int(g.os_event.sum()):2} median OS={k.median_survival_time_:.0f} d")
    cc = conf.assign(k27=(conf.meth_group == "DMG H3 K27-altered").astype(int))
    c = CoxPHFitter().fit(cc[["k27","os_days","os_event"]], "os_days", "os_event")
    print(f"  log-rank p={ps(lr.p_value)}   HR={float(np.exp(c.params_['k27'])):.2f}, p={ps(float(c.summary.loc['k27','p']))}")
    print("\n  CC involvement by molecular subgroup (Fisher):")
    for k_, lab in [("inv_WT","CC+WT"),("inv_TC","CC+TC"),("cross_TC","TC crossing"),
                    ("inv_splenium_WT","splenium WT"),("inv_genu_WT","genu WT")]:
        a, b = int(g1[k_].sum()), int(g0[k_].sum())
        _, p = sps.fisher_exact([[a, len(g1)-a],[b, len(g0)-b]])
        print(f"    {lab:14} K27 {a}/{len(g1)} ({100*a/len(g1):.0f}%)   other {b}/{len(g0)} ({100*b/len(g0):.0f}%)   p={ps(p)}")
else:
    print(f"  not testable: K27 n={len(g1)}, other n={len(g0)}")

# ---------- FIGURES ----------
PED, GBMc = "#c1121f", "#2e86ab"
R = res["OVERALL SURVIVAL"]
fig = plt.figure(figsize=(16.5, 9.2))
gs = fig.add_gridspec(2, 3, width_ratios=[1,1,1.15], hspace=.34, wspace=.26)
for i, (k, lab) in enumerate([("inv_WT","CC involvement (WT)"),("inv_TC","CC involvement (TC)"),
                              ("cross_TC","Midline crossing (TC)"),("inv_splenium_WT","Splenium involvement (WT)")]):
    ax = fig.add_subplot(gs[i//2, i%2])
    d = D.dropna(subset=["os_days"])
    for val, col, nmv in [(1, PED, "involved"), (0, GBMc, "spared")]:
        g = d[d[k] == val]
        if len(g) < 2: continue
        kmf = KaplanMeierFitter(); kmf.fit(g.os_days, g.os_event, label=f"{nmv} (n={len(g)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=col, lw=2.1)
    rr = R[R.label == lab]
    ttl = f"{lab}\nlog-rank {ps(rr.logrank_p.iloc[0])}, HR {rr.hr.iloc[0]:.2f} ({rr.lo.iloc[0]:.2f}-{rr.hi.iloc[0]:.2f})" if len(rr) else lab
    ax.set_title(ttl, fontsize=10)
    ax.set_xlabel("days from diagnosis"); ax.set_ylabel("overall survival probability")
    ax.set_ylim(0, 1.02); ax.grid(alpha=.25); ax.legend(fontsize=8)
axf = fig.add_subplot(gs[:, 2])
y = np.arange(len(R))[::-1]
for i, r in R.iterrows():
    col = PED if r.p < .05 else "#9aa0a6"
    axf.plot([r.lo, r.hi], [y[i]]*2, color=col, lw=2.4); axf.plot(r.hr, y[i], "o", color=col, ms=7.5)
    axf.text(.02, y[i]+.26, r.label, fontsize=8.4, transform=axf.get_yaxis_transform(), ha="left")
    axf.text(.98, y[i]+.26, f"HR {r.hr:.2f}, p={ps(r.p)}", fontsize=7.8,
             transform=axf.get_yaxis_transform(), ha="right", color="#444")
axf.axvline(1, color="black", ls=":", lw=1.3); axf.set_yticks([]); axf.set_xscale("log")
axf.set_xlabel("hazard ratio for death (95% CI)"); axf.grid(alpha=.25, axis="x")
axf.set_ylim(-.7, len(R)-.1); axf.set_title("Unadjusted hazard ratios, overall survival", fontsize=10)
fig.suptitle(f"Exploratory survival analysis, pediatric high-grade glioma "
             f"(n={len(D)} of 114 with outcome data; {int(D.os_event.sum())} deaths)\n"
             f"corpus callosum measured on the JHU atlas; corrected label convention", fontsize=12.5, y=.98)
fig.savefig(OUT/"FigureS3_survival.png", dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
print(f"\nwrote {OUT/'FigureS3_survival.png'}")

# molecular figure
fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
ax = axes[0]
for nm2, g, col in [("DMG H3 K27-altered", g1, PED), ("other pHGG", g0, GBMc)]:
    if len(g) < 2: continue
    k = KaplanMeierFitter(); k.fit(g.os_days, g.os_event, label=f"{nm2} (n={len(g)})")
    k.plot_survival_function(ax=ax, ci_show=True, color=col, lw=2.1)
ax.set_title("Confident methylation calls (score >= 0.9)", fontsize=10)
ax.set_xlabel("days"); ax.set_ylabel("overall survival"); ax.set_ylim(0,1.02); ax.grid(alpha=.25); ax.legend(fontsize=8)
ax = axes[1]
grps = [g for g in D.meth_group.unique() if len(D[D.meth_group == g]) >= 3]
cols = plt.cm.tab10(np.linspace(0, 1, max(len(grps),1)))
for gi, gname in enumerate(grps):
    s = D[D.meth_group == gname]
    k = KaplanMeierFitter(); k.fit(s.os_days, s.os_event, label=f"{gname} (n={len(s)})")
    k.plot_survival_function(ax=ax, ci_show=False, color=cols[gi], lw=2.1)
sub = D[D.meth_group.isin(grps)]
mlr = multivariate_logrank_test(sub.os_days, sub.meth_group, sub.os_event)
ax.set_title(f"All methylation groups\noverall log-rank p={ps(mlr.p_value)}", fontsize=10)
ax.set_xlabel("days"); ax.set_ylabel("overall survival"); ax.set_ylim(0,1.02); ax.grid(alpha=.25); ax.legend(fontsize=7.2)
ax = axes[2]
labs = ["CC+WT","CC+TC","TC cross","genu WT","splen WT"]
keys = ["inv_WT","inv_TC","cross_TC","inv_genu_WT","inv_splenium_WT"]
kv = [100*g1[k].mean() if len(g1) else 0 for k in keys]
ov = [100*g0[k].mean() if len(g0) else 0 for k in keys]
x = np.arange(len(labs)); w = .36
ax.bar(x-w/2, kv, w, color=PED, label=f"K27 (n={len(g1)})")
ax.bar(x+w/2, ov, w, color=GBMc, label=f"other pHGG (n={len(g0)})")
for i in range(len(labs)):
    ax.text(x[i]-w/2, kv[i]+1.5, f"{kv[i]:.0f}%", ha="center", fontsize=8)
    ax.text(x[i]+w/2, ov[i]+1.5, f"{ov[i]:.0f}%", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=9); ax.set_ylim(0, 118)
ax.set_ylabel("% of subgroup"); ax.set_title("Callosal involvement by molecular subgroup", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=.25, axis="y")
fig.suptitle("Molecular subgroup, survival, and callosal involvement", fontsize=12.5, y=1.0)
fig.tight_layout()
fig.savefig(OUT/"FigureS4_molecular.png", dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
print(f"wrote {OUT/'FigureS4_molecular.png'}")
print(f"wrote cc_survival_v3.csv and DataS4_survival_per_patient.csv")
