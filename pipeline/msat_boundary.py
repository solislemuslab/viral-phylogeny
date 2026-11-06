"""
python pipeline/msat_boundary.py X [X ...]

MSA Transformer boundary: the largest N whose random sample (random.seed(42); random.sample(files, N), as in
pipeline/downsample_loci.py) concatenates to at most 1023 columns (1024 positions minus the BOS token).
Prints the column count for N = 1 .. first N above the limit.
"""
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LIMIT = 1023


def locus_length(f):
    with open(f) as h:
        h.readline()
        n = 0
        for line in h:
            if line.startswith(">"):
                break
            n += len(line.strip())
    return n


for X in sys.argv[1:]:
    files = sorted((ROOT / f"{X}45_cleaned").glob("*.fasta"))
    lengths = {f: locus_length(f) for f in files}
    cols = []
    for n in range(1, len(files) + 1):
        random.seed(42)
        cols.append(sum(lengths[f] for f in random.sample(files, n)))
        if cols[-1] > LIMIT and n > 3:
            break
    best = max(n for n, c in enumerate(cols, 1) if c <= LIMIT)
    print(f"{X}: MSA-T boundary N = {best} ({cols[best - 1]} columns); N = {best + 1} has {cols[best]} columns | "
          + " ".join(f"{n}:{c}" for n, c in enumerate(cols, 1)))
