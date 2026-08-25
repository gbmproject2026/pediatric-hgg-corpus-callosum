# Data required, and where to get it

None of the primary data is redistributed here. All of it is publicly available
but licensed to the original providers.

## 1. Pediatric imaging, BraTS-PEDs

The Cancer Imaging Archive, BraTS-PEDs collection.
https://www.cancerimagingarchive.net/collection/brats-peds/

Extract so each case is a folder containing `*-t1c.nii.gz` and `*-seg.nii.gz`:

    <TRAIN_DIR>/BraTS-PED-XXXXX-000/BraTS-PED-XXXXX-000-t1c.nii.gz
    <TRAIN_DIR>/BraTS-PED-XXXXX-000/BraTS-PED-XXXXX-000-seg.nii.gz

The DFCI-BCH-BWH-PEDs-HGG subset has its own TCIA DOI: 10.7937/v8h6-bg25

### Segmentation label convention, pediatric
Four classes (Kazerooni et al., arXiv:2404.15009):

| Value | Region |
|---|---|
| 1 | Enhancing tumor (ET) |
| 2 | Nonenhancing tumor (NET) |
| 3 | Cystic component (CC) |
| 4 | Peritumoral edema (ED) |

This differs from the three-class participant release (NC=1, ED=2, ET=3) and
from the BraTS 2017-2020 adult convention (NCR=1, ED=2, NET=3, ET=4). Neither of
those applies to the four-class files used here.

## 2. Adult imaging, BraTS 2023 adult glioma (GLI)

Synapse, BraTS 2023 challenge data. Requires a Synapse account and auth token.

### Segmentation label convention, adult
Three classes:

| Value | Region |
|---|---|
| 1 | Necrotic core (NCR) |
| 2 | Peritumoral edema (ED) |
| 3 | Enhancing tumor (ET) |

## 3. Clinical and molecular metadata

- `BraTS-PEDs_metadata.tsv`, distributed with BraTS-PEDs
- `cbtn_histologies.tsv`, Children's Brain Tumor Network / OpenPBTA histologies

Both contain patient-level clinical fields and are deliberately excluded from
this repository. Obtain them from the source and place them at the paths set at
the top of the scripts.

## 4. Atlas, included

`atlas_jhu/` holds the JHU ICBM-DTI-81 corpus callosum labels used for all
measurement, already on the MNI152 grid used by the registration target. See
`atlas_jhu/PROVENANCE.md`.
