from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
import json

MSA_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45_cleaned")
OUT_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45_concat")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT = OUT_DIR / "OB45_supermatrix.fasta"

LOG_DIR = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/OB45/logs")
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

with open(LOG_DIR / "taxon_presence_counts.tsv", "w") as out:
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

with open(OUT_DIR / "OB45_partitions.txt", "w") as f:
    f.write("\n".join(partition_lines) + "\n")

summary = {
    "n_alignments": len(files),
    "n_taxa_union": len(taxa_all),
    "supermatrix_length": len(records[0].seq) if records else 0,
    "output": str(OUT),
}
with open(LOG_DIR / "concat_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(summary)
print("Wrote concat/OB45_supermatrix.fasta")
print("Wrote concat/OB45_partitions.txt")
print("Wrote logs/taxon_presence_counts.tsv")
