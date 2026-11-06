"""
Find the largest number of loci Phyloformer can still analyse (in memory) for one dataset.

    python find_phyloformer_boundary.py EB [lo] [hi]

lo is a number of loci known to work and hi one known to fail (default 1 and all loci).
Loci are a random sample of N (random.seed(42); random.sample(files, N), exactly as in the random
downsample_loci_<N>.py scripts) and are concatenated as in concatenate_supermatrix_<N>.py.  Phyloformer is run with the README
commands (infer_alns.py with models/pf.ckpt, then FastME --nni --spr); a probe succeeds when the
tree file is written.  Large probes need several hundred GB of RAM, so probes run one at a time
(binary search) and hold Phyloformer/.big_run.lock so no other big Phyloformer run overlaps.
Every probe (with its peak memory) is logged to <X>45/logs/phyloformer_boundary/probes.tsv.
"""
import random
import subprocess
import sys
import time
from pathlib import Path

from Bio import SeqIO

ROOT = Path("/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia")
PF = ROOT / "Phyloformer"
PYTHON = "/mnt/home/jgao/.conda/envs/phylo/bin/python"
SEED = 42
# Runs a command and prints the peak RSS (kB) of its process tree on the last line
PEAK = ("import resource, subprocess, sys; r = subprocess.run(sys.argv[1:]).returncode; "
        "print('PEAK_KB', resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss); sys.exit(r)")

X = sys.argv[1]
SRC_DIR = ROOT / f"{X}45_cleaned"
WORK = ROOT / f"{X}45/logs/phyloformer_boundary"
WORK.mkdir(parents=True, exist_ok=True)

files = sorted(SRC_DIR.glob("*.fasta"))
lengths = {f: len(next(SeqIO.parse(f, "fasta")).seq) for f in files}
seqs = {f: {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")} for f in files}
taxa = sorted(set().union(*(s.keys() for s in seqs.values())))
lo = int(sys.argv[2]) if len(sys.argv) > 2 else 1
hi = int(sys.argv[3]) if len(sys.argv) > 3 else len(files)


def probe(n):
    random.seed(SEED)
    chosen = sorted(random.sample(files, n))
    d = WORK / f"N{n}"
    (d / "msa").mkdir(parents=True, exist_ok=True)
    (d / "matrix").mkdir(exist_ok=True)
    name = f"{X}45_supermatrix_{n}"
    with open(d / "msa" / f"{name}.fasta", "w") as out:
        for t in taxa:
            out.write(f">{t}\n" + "".join(seqs[f].get(t, "-" * lengths[f]) for f in chosen) + "\n")
    t0 = time.time()
    cmd = (f"flock {PF}/.big_run.lock {PYTHON} -c \"{PEAK}\" {PYTHON} infer_alns.py -o {d}/matrix models/pf.ckpt {d}/msa"
           f" && ./bin/bin_linux/fastme -i {d}/matrix/{name}.phy -o {d}/{name}.nwk --nni --spr")
    with open(d / "phyloformer.out", "w") as log:
        subprocess.run(cmd, shell=True, cwd=PF, stdout=log, stderr=subprocess.STDOUT)
    out = open(d / "phyloformer.out").read()
    ok = (d / f"{name}.nwk").exists()
    peak = next((int(l.split()[1]) for l in out.splitlines() if l.startswith("PEAK_KB")), 0)
    err = next((l.strip()[:120] for l in out.splitlines() if "Error" in l), "")
    (d / "msa" / f"{name}.fasta").unlink()  # large and can be rebuilt
    return ok, sum(lengths[f] for f in chosen), time.time() - t0, peak / 1024 ** 2, err


log = open(WORK / "probes.tsv", "a")
log.write("n_loci\tcolumns\tsuccess\tseconds\tpeak_GB\terror\n")
while hi - lo > 1:
    n = (lo + hi) // 2
    ok, cols, sec, peak, err = probe(n)
    log.write(f"{n}\t{cols}\t{ok}\t{sec:.0f}\t{peak:.0f}\t{err}\n")
    log.flush()
    print(f"N={n} cols={cols} success={ok} {sec:.0f}s peak={peak:.0f}GB {err}", flush=True)
    lo, hi = (n, hi) if ok else (lo, n)

with open(WORK / "boundary.txt", "w") as out:
    random.seed(SEED)
    out.write(f"largest success: {lo} loci ({sum(lengths[f] for f in random.sample(files, lo))} columns)\n")
    out.write(f"smallest failure above it: {hi} loci\n")
print(open(WORK / "boundary.txt").read())
