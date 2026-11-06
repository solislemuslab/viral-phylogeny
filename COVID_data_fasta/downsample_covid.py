"""
Downsample the aligned COVID 300 (+4 outgroups) to N taxa for methods that cannot handle all 304.

    python downsample_covid.py N

Keeps the 4 outgroups and a random sample of N - 4 SARS-CoV-2 sequences (random.seed(42); random.sample),
in their original order, then drops columns that are gaps in every kept sequence.
Writes covid_<N>_aligned.fasta; covid_300_aligned.fasta is not changed.
"""
import random
import sys

from Bio import SeqIO

N = int(sys.argv[1])
OUTGROUPS = {"EPI_ISL_402131", "EPI_ISL_410543", "EPI_ISL_410544", "EPI_ISL_410721"}

recs = list(SeqIO.parse("covid_300_aligned.fasta", "fasta"))
acc = lambda r: r.id.split("|")[1]
ingroup = [r for r in recs if acc(r) not in OUTGROUPS]
assert len(recs) - len(ingroup) == len(OUTGROUPS)

random.seed(42)
chosen = {r.id for r in random.sample(ingroup, N - len(OUTGROUPS))}
kept = [r for r in recs if r.id in chosen or acc(r) in OUTGROUPS]

seqs = [str(r.seq) for r in kept]
cols = [i for i in range(len(seqs[0])) if any(s[i] != "-" for s in seqs)]
with open(f"covid_{N}_aligned.fasta", "w") as out:
    for r, s in zip(kept, seqs):
        out.write(f">{r.id}\n" + "".join(s[i] for i in cols) + "\n")
print(f"covid_{N}_aligned.fasta: {len(kept)} taxa x {len(cols)} columns")
