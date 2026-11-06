from Bio import SeqIO
from Bio.Seq import Seq

INPUT_FASTA = "covid_300_stratified.fasta"
OUTPUT_FASTA = "covid_spike_aa.fasta"

# SARS-CoV-2 Spike coordinates on NC_045512.2
# 1-based inclusive: 21563..25384
SPIKE_START = 21563
SPIKE_END = 25384

records_out = []
kept = 0
skipped = 0

for rec in SeqIO.parse(INPUT_FASTA, "fasta"):
    seq = str(rec.seq).upper().replace("-", "").replace(" ", "")
    
    # skip short genomes
    if len(seq) < SPIKE_END:
        skipped += 1
        continue

    spike_nt = seq[SPIKE_START - 1 : SPIKE_END]

    # length should be divisible by 3
    if len(spike_nt) % 3 != 0:
        skipped += 1
        continue

    aa = str(Seq(spike_nt).translate(to_stop=False))

    # skip sequences with internal stop codons
    # allow terminal stop if present
    if "*" in aa[:-1]:
        skipped += 1
        continue

    aa = aa.rstrip("*")

    rec.seq = Seq(aa)
    records_out.append(rec)
    kept += 1

SeqIO.write(records_out, OUTPUT_FASTA, "fasta")

print(f"Saved {kept} translated Spike amino-acid sequences to {OUTPUT_FASTA}")
print(f"Skipped {skipped} sequences")
