import shutil
from pathlib import Path
from Bio import SeqIO

SRC_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45_cleaned")
OUT_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45_cleaned_short_15")
LOG_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45/logs")
N_KEEP = 15

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(SRC_DIR.glob("*.fasta"))
print(f"Found {len(files)} locus files")

# Keep the N shortest loci (ties broken by file name) so the supermatrix stays under the
# MSA Transformer limit of 1023 columns (1024 positions minus the BOS token).
lengths = {f: len(next(SeqIO.parse(f, "fasta")).seq) for f in files}
chosen = sorted(sorted(files, key=lambda f: (lengths[f], f.name))[:N_KEEP])

for f in chosen:
    shutil.copy(f, OUT_DIR / f.name)

with open(LOG_DIR / "downsample_short_15_selected.txt", "w") as out:
    out.write("\n".join(f"{f.name}\t{lengths[f]}" for f in chosen) + "\n")

print(f"Copied {len(chosen)} files to {OUT_DIR}, total length {sum(lengths[f] for f in chosen)}")
