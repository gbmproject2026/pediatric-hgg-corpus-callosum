# Corpus Callosum Involvement in Pediatric HGG vs Adult GBM

Butterfly glioma, tumor spreading across the corpus callosum, is a familiar sight
in adult glioblastoma. Nobody had measured whether pediatric high-grade gliomas
do the same thing. This project registers 114 pediatric HGG and 1,213 adult
IDH-wildtype GBM cases to MNI space and measures exactly how much tumor sits in
the corpus callosum and where.

The short version of what we found: pediatric tumors reach the corpus callosum
*less* often overall than adult tumors do, but when they get there they push
solid tumor core across the midline more than twice as often (45.0% vs 17.7%,
OR 3.80, P < .001). And the two cohorts favor opposite ends of the structure.
Pediatric involvement sits posteriorly, in the body and splenium. Adult
involvement is anterior, concentrated in the genu, which is what the butterfly
glioma literature has described for years.

## Method

For each case the pipeline registers the post-contrast T1 to the ANTs MNI152
template using SyNRA, with a cost-function mask that excludes the tumor so the
lesion does not drag the alignment. It then warps the expert segmentation into
template space with nearest-neighbor interpolation, which keeps the label values
exact. Overlap with the corpus callosum and its three sub-regions is computed
afterwards, as a separate step.

Keeping those stages separate is deliberate. Registration costs hours of compute
while measurement costs minutes, so any change to the atlas or the label mapping
can be applied by re-running measurement alone.

## Segmentation label conventions

BraTS-PEDs and BraTS-GLI label their segmentations differently, and the values do
not correspond. Reading pediatric label 3 as enhancing tumor, as the adult
convention would suggest, measures the cystic component instead.

| | BraTS-PEDs | BraTS-GLI (adult) |
|---|---|---|
| 1 | Enhancing tumor | Necrotic core |
| 2 | Non-enhancing tumor | Peritumoral edema |
| 3 | Cystic component | Enhancing tumor |
| 4 | Peritumoral edema | not used |

Two further conventions exist in the BraTS family and apply to neither dataset
here: the three-class pediatric participant release (NC=1, ED=2, ET=3) and the
BraTS 2017-2020 adult convention (NCR=1, ED=2, NET=3, ET=4). `docs/METHODS_NOTES.md`
explains how to confirm label identity empirically if you are ever unsure.

Tumor core is taken as all labeled tumor except peritumoral edema in both
cohorts, so labels 1+2+3 for pediatric and 1+3 for adult.

## Setup

Python 3.11 is the easiest path, mostly because antspyx wheels are reliable there.

```bash
python3.11 -m venv venv_ants
source venv_ants/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

If antspyx starts compiling from source, force a wheel with
`pip install antspyx==0.6.3 --only-binary=:all:`.

No script contains a hardcoded path. Everything resolves through
`scripts/config.py`: repository paths come from wherever you cloned it, and the
inputs we cannot redistribute come from environment variables.

```bash
export BRATS_PEDS_TRAIN=/path/to/BraTS-PEDs-v1/Training
export BRATS_GLI_TRAIN=/path/to/BraTS-GLI
export BRATS_PEDS_META=/path/to/BraTS-PEDs_metadata.tsv
export CBTN_HISTOLOGIES=/path/to/histologies.tsv
export CC_WORK=/somewhere/with/space   # optional, defaults to the repo
```

Run `python scripts/config.py` to see what resolved. Anything reported as MISSING
needs setting before you go further. If you would rather not use environment
variables, drop the files under `data/` using the names that command prints.

## Data

Pediatric imaging comes from the BraTS-PEDs collection on TCIA,
https://www.cancerimagingarchive.net/collection/brats-peds/. Extract it so each
case is a folder holding `<id>-t1c.nii.gz` and `<id>-seg.nii.gz`. Adult imaging
comes from the BraTS 2023 challenge on Synapse and needs a free account and an
access token.

Cohort selection also needs `BraTS-PEDs_metadata.tsv` and the CBTN
`histologies.tsv`. Both carry patient-level clinical fields, so neither is
included here. `docs/DATA.md` has the details.

## Measurements

Cases with any overlap in each region. WT is whole tumor, TC is tumor core, ET is
enhancing tumor. TC and ET rows use the enhancing pediatric subgroup (n=80),
whole-tumor rows use all 114.

| Metric | pHGG | Adult GBM (n=1213) |
|--------|------|--------------------|
| CC whole, WT | 81 (71.1%) | 1020 (84.1%) |
| CC whole, TC | 54 (67.5%) | 694 (57.2%) |
| CC whole, ET | 42 (52.5%) | 674 (55.6%) |
| CC genu, WT | 22 (19.3%) | 426 (35.1%) |
| CC genu, TC | 16 (20.0%) | 244 (20.1%) |
| CC body, WT | 74 (64.9%) | 845 (69.7%) |
| CC body, TC | 49 (61.2%) | 444 (36.6%) |
| CC splenium, WT | 61 (53.5%) | 745 (61.4%) |
| CC splenium, TC | 39 (48.8%) | 408 (33.6%) |
| Midline crossing, WT | 50 (43.9%) | 427 (35.2%) |
| Midline crossing, TC | 36 (45.0%) | 215 (17.7%) |
| Midline crossing, ET | 17 (21.2%) | 212 (17.5%) |

Effect sizes and p-values are in `results/cc_results_v3.csv`.

## Registration to MNI152 space

These are real cases warped into template space. Red is the tumor, the contours
are the corpus callosum sub-regions: pink genu, green body, blue splenium. Each
row shows sagittal, coronal and axial views through the callosal centroid. We
picked cases that span the range, including two with no callosal involvement at
all, so the computed flags can be checked against what you can see.

| Case | Cohort | CC whole (WT) | Genu (WT) | Splenium (WT) | TC crossing |
|------|--------|---------------|-----------|---------------|-------------|
| BraTS-PED-00001-000 | pHGG | 18660 | 0 | 11621 | yes |
| BraTS-PED-00138-000 | pHGG | 8293 | 5536 | 0 | yes |
| BraTS-PED-00154-000 | pHGG | 5246 | 0 | 5116 | no |
| BraTS-PED-00136-000 | pHGG | 0 | 0 | 0 | no |
| BraTS-GLI-01432-000 | GBM | 18373 | 8844 | 0 | yes |
| BraTS-GLI-01507-000 | GBM | 0 | 0 | 0 | no |

![QC BraTS-PED-00001-000](qc-samples/QC_BraTS-PED-00001-000.png)
![QC BraTS-PED-00138-000](qc-samples/QC_BraTS-PED-00138-000.png)
![QC BraTS-PED-00154-000](qc-samples/QC_BraTS-PED-00154-000.png)
![QC BraTS-PED-00136-000](qc-samples/QC_BraTS-PED-00136-000.png)
![QC BraTS-GLI-01432-000](qc-samples/QC_BraTS-GLI-01432-000.png)
![QC BraTS-GLI-01507-000](qc-samples/QC_BraTS-GLI-01507-000.png)

The first pediatric case is the pattern in miniature: splenium and body full of
tumor, genu completely clear. The adult case is the mirror image.

## Reproducibility

We registered the pediatric cohort three more times with different random seeds.
The derived classifications agreed 95% to 99% of the time, with Cohen kappa
between 0.90 and 0.95, and prevalence never shifted by more than 3.5 percentage
points. Median whole-tumor Dice against the reference run was 0.90.

Dice is worth reading carefully here. It depends heavily on lesion size, 0.94 in
the largest quartile of tumors against 0.82 in the smallest, so a low value on a
small tumor says more about the metric than about the registration. That is why
we report classification agreement as the headline number. One case did fail
registration outright. Details are in `results/TableS4` through `TableS6` and
`figures/supplemental/FigureS6_reproducibility.png`.

## References

The corpus callosum masks come from the JHU ICBM-DTI-81 white matter atlas,
labels 3, 4 and 5, already on the MNI152 grid we register to. Mori et al., *MRI
Atlas of Human White Matter*, Elsevier 2005, and Hua et al., *NeuroImage*
2008;39(1):336-347. `atlas_jhu/PROVENANCE.md` records the source and the checks
we ran on it.

Pediatric imaging and its label convention are described in Kazerooni et al.,
arXiv:2305.17033 and arXiv:2404.15009, with the DFCI-BCH-BWH-PEDs-HGG subset
under TCIA doi:10.7937/v8h6-bg25. Adult imaging is the BraTS 2023 adult glioma
challenge set. Histologies come from the Children's Brain Tumor Network.
