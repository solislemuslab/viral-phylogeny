from pathlib import Path
from code.extracting import Extractor
from code.build_plm_tree import PlmTree

INPUT_DIR = Path("data/OB45")

OUTPUT_DIR = Path("./OB45_embedding_trees_attention")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

msa_type = "default"
msa_files = sorted(INPUT_DIR.glob("*.fasta"))

print("=" * 80)
print(f"Input directory:  {INPUT_DIR.resolve()}")
print(f"Output directory: {OUTPUT_DIR.resolve()}")
print(f"Found {len(msa_files)} MSA files")
print("=" * 80)

# Run MsaPhylo attention-tree inference for all OB45 files.
for i, msa_file in enumerate(msa_files, start=1):
    # Use the file name without extension as the protein family name.
    # Example: OB_PF00004_1.fasta -> OB_PF00004_1
    prot_family = msa_file.stem

    print("\n" + "=" * 80)
    print(f"[{i}/{len(msa_files)}] Running {prot_family}")
    print(f"MSA file: {msa_file}")
    print("=" * 80)

    try:
        # Step 1: Extract and save column attentions.
        ex = Extractor(prot_family, msa_type)
        ex.get_col_attention()

        # Step 2: Build NJ trees from every attention head.
        plmtree = PlmTree(prot_family, msa_type)
        plmtree.build_attention_tree()

        print(f"Finished {prot_family}")

    except Exception as e:
        print(f"Failed {prot_family}")
        print(f"Error: {e}")
        continue

print("\nAll finished.")
