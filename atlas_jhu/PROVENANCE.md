# atlas_jhu: corpus callosum masks derived from the JHU ICBM-DTI-81 atlas

## Purpose

Corpus callosum masks for measuring tumor overlap, derived from a published
white matter atlas and verified against the registration target.

## Source

JHU ICBM-DTI-81 white matter labels, 48 regions, obtained from the neuroparc
repository:
https://raw.githubusercontent.com/neurodata/neuroparc/master/atlases/label/Human/JHU_space-MNI152NLin6_res-1x1x1.nii.gz

Original atlas reference:
Mori S, Wakana S, van Zijl PCM, Nagae-Poetscher LM. MRI Atlas of Human White
Matter. Elsevier; 2005.
Hua K, Zhang J, Wakana S, et al. Tract probability maps in stereotaxic spaces:
analyses of white matter anatomy and tract-specific quantification.
NeuroImage. 2008;39(1):336-347.

## Labels used

| File | JHU label | Voxels | World y extent |
|---|---|---|---|
| CC_genu.nii.gz | 3, Genu of corpus callosum | 8,851 | +15 to +40 mm |
| CC_body.nii.gz | 4, Body of corpus callosum | 13,711 | -30 to +19 mm |
| CC_splenium.nii.gz | 5, Splenium of corpus callosum | 12,729 | -68 to -31 mm |
| CC_mask.nii.gz | 3 or 4 or 5 | 35,291 | -68 to +40 mm |

Anterior-to-posterior ordering (genu +26 mm, body -5 mm, splenium -42 mm)
matches published MNI coordinates for these structures.

## Grid

The atlas is already 182 x 218 x 182 with affine
[[-1,0,0,90],[0,1,0,-126],[0,0,1,-72]], identical to the ANTs MNI template
returned by `ants.get_ants_data('mni')` and therefore identical to the grid of
every registered segmentation in `cc_results/registrations/` and
`cc_results_gbm/registrations/`. No resampling or interpolation is required.

Mean template T1 intensity inside the new CC mask is 1.19x whole-brain mean,
consistent with white matter.

## Verification performed

- Grid and affine equality with the registration target: confirmed.
- Subregion anterior-posterior ordering against published MNI coordinates: confirmed.
- Bounding-box fill ratios well below 1.00, i.e. shaped structures: confirmed.
- Template intensity ratio above 1.0, i.e. white matter: confirmed.
