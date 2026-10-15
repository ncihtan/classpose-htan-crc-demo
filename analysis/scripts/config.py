"""Single source of truth for analysis paths. Every value is env-var overridable."""
import os
from pathlib import Path

def _p(var: str, default: Path) -> Path:
    return Path(os.environ.get(var, default)).expanduser().resolve()

ROOT            = _p("CLASSPOSE_ROOT", Path(__file__).resolve().parents[1])   # analysis/
SAMPLE_DATA_DIR = _p("SAMPLE_DATA_DIR", ROOT / "sample_data")   # zarrs + tifs + mIF csv (Synapse/S3, gitignored)
METADATA_DIR    = _p("METADATA_DIR", ROOT / "metadata")
DATA_DIR        = _p("DATA_DIR", ROOT / "data")
FIGURES_DIR     = _p("FIGURES_DIR", ROOT / "figures")
MIF_CSV         = _p("MIF_CSV", SAMPLE_DATA_DIR / "P37_S29_A24_C59kX_E15__at__20220106_014304_946511_cellRingMask.csv")
S3_BASE         = os.environ.get("CLASSPOSE_S3_BASE",
                                 "s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model")
AWS_PROFILE     = os.environ.get("AWS_PROFILE")   # None -> default credential chain
