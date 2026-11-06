from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import json
import sys

# python concatenate_supermatrix_n.py N  ->  EB45_concat_N/EB45_supermatrix_N.fasta
N = int(sys.argv[1])

MSA_DIR = Path(f"/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/EB45_cleaned_{N}")
OUT_DIR = Path(f"/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/EB45_concat_{N}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT = OUT_DIR / f"EB45_supermatrix_{N}.fasta"

LOG_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/EB45/logs")
LOG_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(MSA_DIR.glob("*.fasta"))

# First pass: collect taxa and alignment lengths
taxa_all = set()
aln_info = []

for f in files:
    records = list(SeqIO.parse(f, "fasta"))
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        raise ValueError(f"{f} has inconsistent sequence lengths: {lengths}")

    L = list(lengths)[0]
    taxa = [r.id for r in records]
    taxa_all.update(taxa)
    aln_info.append((f, L, set(taxa)))

print(f"Number of alignments: {len(files)}")
print(f"Number of unique taxa across all alignments: {len(taxa_all)}")

# Check how often each taxon appears
taxon_counts = Counter()
for _, _, taxa in aln_info:
    for t in taxa:
        taxon_counts[t] += 1

with open(LOG_DIR / f"taxon_presence_counts_{N}.tsv", "w") as out:
    out.write("taxon\tcount\n")
    for taxon, count in taxon_counts.most_common():
        out.write(f"{taxon}\t{count}\n")

# Use union of taxa, missing OG segments filled with gaps
supermatrix = {taxon: [] for taxon in sorted(taxa_all)}
partition_lines = []
start = 1

for idx, (f, L, taxa) in enumerate(aln_info, start=1):
    recs = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}

    for taxon in supermatrix:
        if taxon in recs:
            supermatrix[taxon].append(recs[taxon])
        else:
            supermatrix[taxon].append("-" * L)

    end = start + L - 1
    partition_lines.append(f"OG{idx} = {start}-{end} # {f.name}")
    start = end + 1

records = [
    SeqRecord(Seq("".join(parts)), id=taxon, description="")
    for taxon, parts in supermatrix.items()
]

SeqIO.write(records, OUT, "fasta")

with open(OUT_DIR / f"EB45_partitions_{N}.txt", "w") as f:
    f.write("\n".join(partition_lines) + "\n")

summary = {
    "n_alignments": len(files),
    "n_taxa_union": len(taxa_all),
    "supermatrix_length": len(records[0].seq) if records else 0,
    "output": str(OUT),
}
with open(LOG_DIR / f"concat_summary_{N}.json", "w") as f:
    json.dump(summary, f, indent=2)

print(summary)
print(f"Wrote EB45_concat_{N}/EB45_supermatrix_{N}.fasta")
print(f"Wrote EB45_concat_{N}/EB45_partitions_{N}.txt")
print(f"Wrote logs/taxon_presence_counts_{N}.tsv")
