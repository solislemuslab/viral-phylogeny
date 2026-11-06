import random
import shutil
from pathlib import Path

SRC_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45_cleaned")
OUT_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45_cleaned_100")
LOG_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45/logs")
N_KEEP = 100
SEED = 42 

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(SRC_DIR.glob("*.fasta"))
print(f"Found {len(files)} locus files")

random.seed(SEED)
chosen = sorted(random.sample(files, N_KEEP))

for f in chosen:
    shutil.copy(f, OUT_DIR / f.name)

with open(LOG_DIR / "downsample_100_selected.txt", "w") as out:
    out.write("\n".join(f.name for f in chosen) + "\n")

print(f"Copied {len(chosen)} files to {OUT_DIR}")
