from Bio import SeqIO
import pandas as pd
import random

# --------------------------
# SETTINGS
# --------------------------
FASTA_FILE = "covid_unique.fasta"

META_FILE_1 = "gisaid_hcov-19_to_2020_03_30.csv"
META_FILE_2 = "gisaid_hcov-19_after_2020_03_30.csv"

OUTPUT_FASTA = "covid_300_stratified.fasta"
OUTPUT_META = "covid_300_stratified_metadata.csv"

TOTAL_N = 300
random.seed(42)

# --------------------------
# Helper functions
# --------------------------
def load_metadata(path: str) -> pd.DataFrame:
    # Skip the first title row, use the second row as header
    df = pd.read_csv(path, skiprows=1)

    # Fallback in case export is tab-separated or semicolon-separated
    if len(df.columns) == 1:
        df = pd.read_csv(path, skiprows=1, sep="\t")
    if len(df.columns) == 1:
        df = pd.read_csv(path, skiprows=1, sep=";")

    df.columns = [str(c).strip() for c in df.columns]
    return df

def parse_continent_from_location(location: str):
    """
    GISAID 'Location' is usually like:
    'Asia / China / Wuhan'
    'Europe / Italy / Lombardy'
    """
    if pd.isna(location):
        return None
    parts = [x.strip() for x in str(location).split("/")]
    if len(parts) == 0:
        return None
    return parts[0]

def normalize_accession(x):
    if pd.isna(x):
        return None
    return str(x).strip()

# --------------------------
# Load and merge metadata
# --------------------------
meta1 = load_metadata(META_FILE_1)
meta2 = load_metadata(META_FILE_2)

print("meta1 columns:", meta1.columns.tolist())
print("meta2 columns:", meta2.columns.tolist())

meta = pd.concat([meta1, meta2], ignore_index=True)

# --------------------------
# Metadata columns
# --------------------------
ID_COL = "Accession ID"
LOCATION_COL = "Location"

required_cols = [ID_COL, LOCATION_COL]
missing = [c for c in required_cols if c not in meta.columns]
if missing:
    raise ValueError(
        f"Missing columns in metadata: {missing}\n"
        f"Available columns are: {list(meta.columns)}"
    )

meta = meta[[ID_COL, LOCATION_COL]].copy()
meta[ID_COL] = meta[ID_COL].apply(normalize_accession)
meta["continent"] = meta[LOCATION_COL].apply(parse_continent_from_location)

# Drop missing / duplicate metadata rows
meta = meta.dropna(subset=[ID_COL, "continent"])
meta = meta.drop_duplicates(subset=[ID_COL])

# --------------------------
# Read FASTA
# --------------------------
records = list(SeqIO.parse(FASTA_FILE, "fasta"))

# --------------------------
# Match FASTA records to metadata
# Assumes accession ID appears somewhere in FASTA header
# e.g. >...|EPI_ISL_4212|...
# --------------------------
matched = []
unmatched = 0

metadata_ids = set(meta[ID_COL].tolist())

for rec in records:
    header = rec.description
    found_accession = None

    for acc in metadata_ids:
        if acc in header:
            found_accession = acc
            break

    if found_accession is not None:
        continent = meta.loc[meta[ID_COL] == found_accession, "continent"].iloc[0]
        matched.append((rec, found_accession, continent))
    else:
        unmatched += 1

print(f"Matched records: {len(matched)}")
print(f"Unmatched records: {unmatched}")

if len(matched) == 0:
    raise ValueError(
        "No FASTA headers matched metadata accession IDs. "
        "Check whether accession IDs appear in your FASTA headers."
    )

df = pd.DataFrame(matched, columns=["record", "accession", "continent"])

# --------------------------
# Stratified sampling by continent
# --------------------------
counts = df["continent"].value_counts()
proportions = counts / counts.sum()

samples_per_group = (proportions * TOTAL_N).round().astype(int)

# Fix rounding so total exactly = TOTAL_N
diff = TOTAL_N - samples_per_group.sum()
if diff != 0:
    samples_per_group.iloc[0] += diff

selected_records = []
selected_rows = []

for continent, n in samples_per_group.items():
    subset = df[df["continent"] == continent].copy()

    if len(subset) <= n:
        chosen = subset
    else:
        chosen = subset.sample(n=n, random_state=42)

    selected_records.extend(chosen["record"].tolist())
    selected_rows.append(chosen)

selected_df = pd.concat(selected_rows, ignore_index=True)

# If somehow not exactly TOTAL_N due to edge cases, trim or top up
if len(selected_df) > TOTAL_N:
    selected_df = selected_df.sample(n=TOTAL_N, random_state=42)
    selected_records = selected_df["record"].tolist()

elif len(selected_df) < TOTAL_N:
    used_acc = set(selected_df["accession"])
    remaining = df[~df["accession"].isin(used_acc)]
    need = TOTAL_N - len(selected_df)
    if len(remaining) > 0:
        extra = remaining.sample(n=min(need, len(remaining)), random_state=42)
        selected_df = pd.concat([selected_df, extra], ignore_index=True)
        selected_records = selected_df["record"].tolist()

# --------------------------
# Save outputs
# --------------------------
SeqIO.write(selected_records, OUTPUT_FASTA, "fasta")
selected_df[["accession", "continent"]].to_csv(OUTPUT_META, index=False)

print(f"Saved {len(selected_records)} sequences to {OUTPUT_FASTA}")
print(f"Saved metadata to {OUTPUT_META}")
print("\nSelected counts by continent:")
print(selected_df["continent"].value_counts())
