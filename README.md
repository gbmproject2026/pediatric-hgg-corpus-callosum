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

Prevalence ratios with raw and Benjamini-Hochberg adjusted p-values are in
`results/cc_results_v3.csv`.

## Statistics


| Methods | Question |
|---|---|
| Fisher exact with Benjamini-Hochberg | Do the cohorts differ in the proportion involved, across 15 outcomes |
| Prevalence ratio with 95% CI | How large is the difference |
| Mann-Whitney U | Do tumor volumes, and the fraction of corpus callosum involved, differ |
| Modified Poisson regression | Does the difference survive adjustment for tumor volume |

 

Tumor core midline crossing is the pre-specified primary outcome and is reported
unadjusted; the other 14 carry a Benjamini-Hochberg adjusted p-value.

## Results

pHGG n=114 (80 enhancing), adult GBM n=1213. 5/15 outcomes survive
correction:

| Finding | pHGG | GBM | PR (95% CI) | P | P (BH) |
|---|---|---|---|---|---|
| **Midline crossing, TC (>20/side)** | 45.0% | 17.7% | 2.54 (1.94-3.32) | <.001 | <.001 |
| **Body involvement, TC** | 61.2% | 36.6% | 1.67 (1.38-2.02) | <.001 | <.001 |
| **Genu involvement, WT** | 19.3% | 35.1% | 0.56 (0.38-0.81) | <.001 | .003 |
| **CC involvement, WT** | 71.1% | 84.1% | 0.84 (0.75-0.95) | <.001 | .004 |
| **Splenium involvement, TC** | 48.8% | 33.6% | 1.45 (1.14-1.84) | .007 | .022 |

 
## Results
Tumor core midline crossing is 2.5-fold more frequent in pediatric high-grade glioma: (36/80, 45.0%) than in adult glioblastoma (215/1213, 17.7%; PR 2.54, 95% CI 1.94 to 3.32; P < .001, BH-adjusted P < .001). The association persists after adjustment for whole-tumor volume (adjusted PR 4.23, 95% CI 3.20 to 5.60) 

Sub-regional involvement shows opposing anterior-posterior gradients: Pediatric tumor core involves the body (61.2% vs 36.6%; PR 1.67, 95% CI 1.38 to 2.02; P < .001) and splenium (48.8% vs 33.6%; PR 1.45, 95% CI 1.14 to 1.84; P = .007) more frequently, whereas genu involvement by whole tumor is less frequent (19.3% vs 35.1%; PR 0.56, 95% CI 0.38 to 0.81; P < .001). All three comparisons survive correction. The anterior predominance in the adult cohort is concordant with existing literature on butterfly glioblastomas.

Callosal burden is greater in pediatric cases among those involved: Analyzed as a continuous fraction, pediatric tumor core occupies 8.7% of the corpus callosum versus 5.5% in adults (Mann-Whitney U, P = .0009), indicating a difference in extent as well as in frequency, and independent of any threshold.

The higher overall involvement in adults is attributable to tumor volume: Adult glioblastomas are 2.45-fold larger (Mann-Whitney U, P < .001), and the whole-tumor difference does not persist after volume adjustment (adjusted PR 1.08, 95% CI 0.97 to 1.20; P = .19).

### Survival

Exploratory, in the 49 patients with outcome data (42 deaths, median OS 527
days). Callosal involvement is not associated with overall survival for any
compartment or sub-region, including tumor core midline crossing (HR 1.17,
95% CI 0.62-2.20, P = .64). Event-free survival shows no association either.

Methylation data available for 37 patients: Survival did not differ between H3 K27-altered tumors and
other pediatric high-grade gliomas (HR 1.33, P = .63, n=17). No K27-altered tumor
involved the genu (0/12) compared with 3 of 5 other pediatric tumors (P = .015),
consistent with the posterior pattern seen in the findings.

![Survival](figures/supplemental/FigureS3_survival.png)

 

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


