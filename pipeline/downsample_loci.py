"""
python pipeline/downsample_loci.py X N  ->  <X>45_cleaned_<N>/ with N randomly selected loci of <X>45_cleaned/

Same sampling as OB45/downsample_loci_7.py: random.seed(42); random.sample(sorted files, N).
The selected file names are written to <X>45/logs/downsample_<N>_selected.txt.
"""
import random
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
X, N_KEEP = sys.argv[1], int(sys.argv[2])
SRC_DIR = ROOT / f"{X}45_cleaned"
OUT_DIR = ROOT / f"{X}45_cleaned_{N_KEEP}"
LOG_DIR = ROOT / f"{X}45/logs"
SEED = 42

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(SRC_DIR.glob("*.fasta"))
print(f"Found {len(files)} locus files")

random.seed(SEED)
chosen = sorted(random.sample(files, N_KEEP))

for f in chosen:
    shutil.copy(f, OUT_DIR / f.name)

with open(LOG_DIR / f"downsample_{N_KEEP}_selected.txt", "w") as out:
    out.write("\n".join(f.name for f in chosen) + "\n")

print(f"Copied {len(chosen)} files to {OUT_DIR}")
