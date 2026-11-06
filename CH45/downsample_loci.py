import random
import shutil
import sys
from pathlib import Path

# python downsample_loci.py N  ->  CH45_cleaned_N with N randomly selected loci (same sampling as OB45/downsample_loci_7.py)
SRC_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/CH45_cleaned")
N_KEEP = int(sys.argv[1])
OUT_DIR = Path(f"/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/CH45_cleaned_{N_KEEP}")
LOG_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/CH45/logs")
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
