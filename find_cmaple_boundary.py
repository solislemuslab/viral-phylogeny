"""
Find the largest number of loci CMAPLE can still analyse for one dataset.

    python find_cmaple_boundary.py AS [lo] [hi]

lo is a number of loci known to work and hi one known to fail (default 0 and all loci).

Loci are a random sample of N (random.seed(42); random.sample(files, N)) exactly as in the random
downsample_loci_<N>.py scripts, and are concatenated in file-name order exactly as in
concatenate_supermatrix_<N>.py, so a probe of N loci is the same supermatrix the pipeline later builds for N.  CMAPLE is run with the README settings
(cmaple-aa, seed 697640); a probe succeeds when CMAPLE writes a tree file.

Each round probes several N in parallel between the largest success (lo) and the smallest failure (hi),
geometrically spaced because the boundary is usually small, until hi == lo + 1.  Success is not monotonic in N (a sample
can fail while a larger sample containing it succeeds), so a success above a failure is possible and is reported.
Every probe is logged to <X>45/logs/cmaple_boundary/probes.tsv and the answer to boundary.txt.
"""
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import SeqIO

ROOT = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia")
CMAPLE = ROOT / "CMAPLE/cmaple-1.1.0-Linux-intel/bin/cmaple-aa"
SEED = "697640"  # CMAPLE seed from the README
SAMPLE_SEED = 42
N_PARALLEL = 6

X = sys.argv[1]
SRC_DIR = ROOT / f"{X}45_cleaned"
WORK = ROOT / f"{X}45/logs/cmaple_boundary"
WORK.mkdir(parents=True, exist_ok=True)

files = sorted(SRC_DIR.glob("*.fasta"))
lengths = {f: len(next(SeqIO.parse(f, "fasta")).seq) for f in files}
seqs = {f: {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")} for f in files}
taxa = sorted(set().union(*(s.keys() for s in seqs.values())))


def probe(n):
    random.seed(SAMPLE_SEED)
    chosen = sorted(random.sample(files, n))
    d = WORK / f"N{n}"
    d.mkdir(exist_ok=True)
    aln = d / f"{X}45_supermatrix_{n}.fasta"
    with open(aln, "w") as out:
        for t in taxa:
            out.write(f">{t}\n" + "".join(seqs[f].get(t, "-" * lengths[f]) for f in chosen) + "\n")
    t0 = time.time()
    with open(d / "cmaple.out", "w") as log:
        subprocess.run([str(CMAPLE), "-aln", aln.name, "--prefix", f"{X}45_supermatrix_{n}", "-seed", SEED],
                       cwd=d, stdout=log, stderr=subprocess.STDOUT)
    ok = (d / f"{X}45_supermatrix_{n}.treefile").exists()
    err = next((l.strip() for l in open(d / "cmaple.out") if l.startswith("ERROR")), "")
    aln.unlink()  # the probe alignments are large and can be rebuilt
    return n, ok, sum(lengths[f] for f in chosen), time.time() - t0, err


done = {}
lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
hi = int(sys.argv[3]) if len(sys.argv) > 3 else len(files)  # the full supermatrix is run by run_<x>45.sh
log = open(WORK / "probes.tsv", "a")
log.write("n_loci\tcolumns\tsuccess\tseconds\terror\n")
while hi - lo > 1:
    a = max(lo, 1)
    pts = sorted({min(hi - 1, max(lo + 1, round(a * (hi / a) ** (i / (N_PARALLEL + 1)))))
                  for i in range(1, N_PARALLEL + 1)} - done.keys())
    with ThreadPoolExecutor(len(pts)) as ex:
        for n, ok, cols, sec, err in ex.map(probe, pts):
            done[n] = ok
            log.write(f"{n}\t{cols}\t{ok}\t{sec:.0f}\t{err}\n")
            log.flush()
            print(f"N={n} cols={cols} success={ok} {sec:.0f}s {err}", flush=True)
    lo = max([n for n, ok in done.items() if ok and n < hi] + [lo])
    hi = min([n for n, ok in done.items() if not ok and n > lo] + [hi])

non_monotone = sorted(n for n, ok in done.items() if ok and n > hi)
with open(WORK / "boundary.txt", "w") as out:
    random.seed(SAMPLE_SEED)
    out.write(f"largest success: {lo} loci ({sum(lengths[f] for f in random.sample(files, lo))} columns)\n")
    out.write(f"smallest failure above it: {hi} loci\n")
    if non_monotone:
        out.write(f"successes above the failure (non-monotone): {non_monotone}\n")
print(open(WORK / "boundary.txt").read())
