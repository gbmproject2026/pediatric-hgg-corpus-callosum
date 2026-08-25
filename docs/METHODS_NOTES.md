# Methodological notes

Technical requirements for anyone running this pipeline.

## Corpus callosum masks

Measurement uses the JHU ICBM-DTI-81 white matter atlas, labels 3 (genu),
4 (body) and 5 (splenium). See `atlas_jhu/PROVENANCE.md`.

`scripts/cc_build_atlas.py` verifies four properties before writing the masks,
and exits if any fails:

- bounding-box fill ratio below 0.9, confirming an anatomically shaped structure
- mean template T1 intensity above whole-brain mean, confirming white matter
- anterior-to-posterior ordering of genu, body and splenium
- the three sub-regions partitioning the whole mask exactly, with no overlap

Sub-region masks should be checked in a sagittal plane. In an axial view the
corpus callosum appears as two separate limbs, which makes anterior and posterior
structures difficult to distinguish visually.

## Segmentation label conventions

BraTS-PEDs, four classes (Kazerooni et al., arXiv:2404.15009):

| Value | Region |
|---|---|
| 1 | Enhancing tumor (ET) |
| 2 | Nonenhancing tumor (NET) |
| 3 | Cystic component (CC) |
| 4 | Peritumoral edema (ED) |

BraTS-GLI adult, three classes:

| Value | Region |
|---|---|
| 1 | Necrotic core (NCR) |
| 2 | Peritumoral edema (ED) |
| 3 | Enhancing tumor (ET) |

Two further conventions exist in the BraTS family and apply to neither dataset
used here: the three-class pediatric participant release (NC=1, ED=2, ET=3), and
the BraTS 2017-2020 adult convention (NCR=1, ED=2, NET=3, ET=4).

Label identity can be confirmed from the imaging alone if a dataset is
undocumented. Enhancing tumor is the class with the greatest median signal
increase from pre-contrast to post-contrast T1, a within-case comparison that
needs no normalization. The cystic component is the class that is simultaneously
hyperintense on T2 and hypointense on post-contrast T1. Applied to the adult
cohort, where the convention is known, this recovers the enhancing label in 97%
of cases.

Tumor core is defined as all labeled tumor except peritumoral edema in both
cohorts: labels 1+2+3 for pediatric and 1+3 for adult.

## Cohort selection

The CBTN histology table records one row per biospecimen, up to 198 rows per
participant, and the primary row is frequently unpopulated. `cc_cohort_v2.py`
scans every row for a participant rather than the first only.

The pontine exclusion is applied on `pathology_diagnosis`. The `primary_site`
field is multi-valued, for example "Thalamus;Brain Stem- Midbrain/Tectum", so
filtering on it would also remove eligible thalamic cases.

## Reproducibility

Repeat registration with independent random seeds gives 95% to 99% agreement on
the derived classifications (Cohen kappa 0.90 to 0.95).

Voxelwise Dice is volume-dependent and is not a registration-quality threshold on
its own: it is substantially lower for small lesions at equal registration
accuracy. Where the outcome is a binary classification, report classification
agreement rather than Dice.
