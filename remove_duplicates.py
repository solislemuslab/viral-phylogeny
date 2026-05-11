from Bio import SeqIO

seen = set()
unique = []

for record in SeqIO.parse("covid_all.fasta", "fasta"):
    seq = str(record.seq)
    if seq not in seen:
        seen.add(seq)
        unique.append(record)

SeqIO.write(unique, "covid_unique.fasta", "fasta")
print("Final unique sequences:", len(unique))
