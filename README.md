# Corpus Callosum Involvement in Pediatric HGG vs Adult GBM



## Method

For each case the pipeline registers the post-contrast T1 to the ANTs MNI152
template using SyNRA, with a cost-function mask that excludes the tumor. It then warps the expert segmentation into
template space with nearest-neighbor interpolation, which keeps the label values
exact. Overlap with the corpus callosum and its three sub-regions is then computed.

## Segmentation label conventions


| | BraTS-PEDs | BraTS-GLI (adult) |
|---|---|---|
| 1 | Enhancing tumor | Necrotic core |
| 2 | Non-enhancing tumor | Peritumoral edema |
| 3 | Cystic component | Enhancing tumor |
| 4 | Peritumoral edema | not used |


## Data

Pediatric imaging comes from the BraTS-PEDs collection on TCIA,
https://www.cancerimagingarchive.net/collection/brats-peds/. Extract it so each
case is a folder holding `<id>-t1c.nii.gz` and `<id>-seg.nii.gz`. Adult imaging
comes from the BraTS 2023 challenge on Synapse and needs a free account and an
access token.

Cohort selection also needs `BraTS-PEDs_metadata.tsv` and the CBTN
`histologies.tsv`.

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



## Reproducibility

We registered the pediatric cohort three more times with different random seeds.
The derived classifications agreed 95% to 99% of the time, with Cohen kappa
between 0.90 and 0.95, and prevalence never shifted by more than 3.5 percentage
points. Median whole-tumor Dice against the reference run was 0.90.


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
