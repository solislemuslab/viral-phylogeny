import shutil
from pathlib import Path
from Bio import SeqIO

SRC_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/EB45_cleaned")
OUT_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/EB45_cleaned_26")
LOG_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/EB45/logs")
N_KEEP = 26

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(SRC_DIR.glob("*.fasta"))
print(f"Found {len(files)} locus files")

# Keep the N shortest loci (ties broken by file name) (the boundary search in find_cmaple_boundary.py
# and the MSA Transformer limit of 1023 columns both use this same ordering).
lengths = {f: len(next(SeqIO.parse(f, "fasta")).seq) for f in files}
chosen = sorted(sorted(files, key=lambda f: (lengths[f], f.name))[:N_KEEP])

for f in chosen:
    shutil.copy(f, OUT_DIR / f.name)

with open(LOG_DIR / "downsample_26_selected.txt", "w") as out:
    out.write("\n".join(f"{f.name}\t{lengths[f]}" for f in chosen) + "\n")

print(f"Copied {len(chosen)} files to {OUT_DIR}, total length {sum(lengths[f] for f in chosen)}")
