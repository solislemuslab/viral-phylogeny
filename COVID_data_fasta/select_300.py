from Bio import SeqIO
import random

records = list(SeqIO.parse("covid_unique.fasta", "fasta"))

random.seed(42)
selected = random.sample(records, 300)

SeqIO.write(selected, "covid_300.fasta", "fasta")

print("Selected 300 sequences")
