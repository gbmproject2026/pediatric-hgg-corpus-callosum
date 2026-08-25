"""
config.py
Single place where every path is resolved. Import this instead of hardcoding.

Repository paths are derived from this file's location, so a clone works with no
edits. The three external inputs are not redistributable and are read from
environment variables, with a `data/` fallback inside the repo.

Environment variables
  BRATS_PEDS_TRAIN   directory of BraTS-PEDs case folders (default data/BraTS-PEDs-v1/Training)
  BRATS_GLI_TRAIN    directory of BraTS-GLI case folders  (default data/BraTS-GLI)
  BRATS_PEDS_META    BraTS-PEDs_metadata.tsv             (default data/BraTS-PEDs_metadata.tsv)
  CBTN_HISTOLOGIES   CBTN histologies.tsv                (default data/histologies.tsv)
  CC_WORK            where warped segmentations are written (default repo root)

Example
  export BRATS_PEDS_TRAIN=/data/BraTS-PEDs-v1/Training
  export BRATS_GLI_TRAIN=/data/BraTS-GLI
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# committed, always present
ATLAS = REPO / "atlas_jhu"
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
FIGURES_SUP = FIGURES / "supplemental"
QC_SAMPLES = REPO / "qc-samples"

# derived, created by the pipeline
WORK = Path(os.environ.get("CC_WORK", REPO))
PEDS_REG = WORK / "cc_results_peds" / "registrations"
GLI_REG = WORK / "cc_results_gbm" / "registrations"


def seed_reg(seed):
    """Warped segmentations from a repeat registration with a given seed."""
    return WORK / f"cc_results_seed{seed}" / "registrations"


# external inputs, not redistributed. See docs/DATA.md
def _ext(var, default):
    return Path(os.environ.get(var, REPO / "data" / default))


PEDS_TRAIN = _ext("BRATS_PEDS_TRAIN", "BraTS-PEDs-v1/Training")
GLI_TRAIN = _ext("BRATS_GLI_TRAIN", "BraTS-GLI")
PEDS_META = _ext("BRATS_PEDS_META", "BraTS-PEDs_metadata.tsv")
CBTN = _ext("CBTN_HISTOLOGIES", "histologies.tsv")

# pipeline artefacts
COHORT = RESULTS / "cc_cohort_v2.txt"
COHORT_REASONS = RESULTS / "cc_cohort_v2_reasons.csv"
GBM_CASES = RESULTS / "gbm_cases.txt"
SIDES = RESULTS / "cc_sides_v3.csv"
OUTCOMES = RESULTS / "cc_results_v3.csv"
SURVIVAL = RESULTS / "cc_survival_v3.csv"

for _d in (RESULTS, FIGURES, FIGURES_SUP, QC_SAMPLES):
    _d.mkdir(parents=True, exist_ok=True)


def require(path, what, hint="see docs/DATA.md"):
    """Fail with an actionable message rather than a stack trace."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"missing {what}: {p}\n  {hint}")
    return p


if __name__ == "__main__":
    print(f"repo         {REPO}")
    print(f"atlas        {ATLAS}          {'ok' if ATLAS.exists() else 'MISSING'}")
    print(f"results      {RESULTS}")
    print(f"work         {WORK}")
    print(f"peds regs    {PEDS_REG}      {'ok' if PEDS_REG.exists() else 'not yet built'}")
    print(f"gli regs     {GLI_REG}       {'ok' if GLI_REG.exists() else 'not yet built'}")
    print()
    for n, p in [("BRATS_PEDS_TRAIN", PEDS_TRAIN), ("BRATS_GLI_TRAIN", GLI_TRAIN),
                 ("BRATS_PEDS_META", PEDS_META), ("CBTN_HISTOLOGIES", CBTN)]:
        print(f"{n:18} {p}   {'ok' if p.exists() else 'MISSING, set the env var'}")
