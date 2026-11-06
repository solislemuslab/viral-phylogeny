"""
python pipeline/preprocess.py X  ->  <X>45_cleaned/: one cleaned FASTA per alignment in <X>45/ (raw .afa files)

Same as OB45/preprocess_ob45.py and AG45/preprocess_ag45.py (taxon code = first word of the header, invalid
characters -> X), plus B/Z -> X because CMAPLE and Phyloformer reject them (only CH45 has any).
"""
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
X = sys.argv[1]
RAW_DIR = ROOT / f"{X}45"
OUT_DIR = ROOT / f"{X}45_cleaned"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_DIR = ROOT / f"{X}45/logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO-*?.")
EXTS = {".fa", ".fas", ".fasta", ".aln", ".faa", ".afa", ".txt"}

def clean_taxon_id(description):
    """
    Example original header:
    >SOLCM H8KWX4.1

    We keep only SOLCM as the taxon/species code.
    This is important for concatenating the alignments.
    """
    first_token = description.strip().split()[0]
    first_token = re.sub(r"[^A-Za-z0-9_.|:-]", "_", first_token)
    return first_token

summary = []
bad_files = []

files = [p for p in RAW_DIR.rglob("*") if p.is_file() and p.suffix.lower() in EXTS]

print(f"Found candidate MSA files: {len(files)}")

for path in sorted(files):
    records = list(SeqIO.parse(path, "fasta"))

    if not records:
        bad_files.append((str(path), "no_records"))
        continue

    cleaned = []
    seen = set()

    for idx, rec in enumerate(records, start=1):
        taxon_id = clean_taxon_id(rec.description)

        if taxon_id in seen:
            taxon_id = f"{taxon_id}_{idx}"
        seen.add(taxon_id)

        seq = str(rec.seq).upper().replace(" ", "").replace("\n", "")
        seq = "".join(ch if ch in VALID_AA else "X" for ch in seq)
        # Some datasets (CH45) have a few ambiguity codes (B, Z) that CMAPLE and Phyloformer reject; treat them as unknown
        seq = seq.replace("B", "X").replace("Z", "X")

        if len(seq) == 0:
            continue

        cleaned.append(SeqRecord(Seq(seq), id=taxon_id, description=""))

    if len(cleaned) < 4:
        bad_files.append((str(path), f"too_few_sequences:{len(cleaned)}"))
        continue

    lengths = {len(r.seq) for r in cleaned}

    if len(lengths) != 1:
        bad_files.append((str(path), f"unaligned_lengths:{sorted(lengths)}"))
        continue

    out_name = path.stem
    out_name = re.sub(r"[^A-Za-z0-9_.-]", "_", out_name)
    out_path = OUT_DIR / f"{out_name}.fasta"

    SeqIO.write(cleaned, out_path, "fasta")

    summary.append({
        "input": str(path),
        "output": str(out_path),
        "nseq": len(cleaned),
        "length": list(lengths)[0],
        "taxa": [r.id for r in cleaned],
    })

with open(LOG_DIR / "preprocess_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

with open(LOG_DIR / "bad_files.txt", "w") as f:
    for item in bad_files:
        f.write("\t".join(item) + "\n")

print(f"Cleaned MSA files: {len(summary)}")
print(f"Bad files: {len(bad_files)}")
print("Wrote logs/preprocess_summary.json and logs/bad_files.txt")
