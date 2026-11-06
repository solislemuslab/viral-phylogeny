"""
MsaPhylo attention trees for one supermatrix (copy of MsaPhylo/run_ag45_attention_trees_18.py with the dataset as arguments).
Run from the MsaPhylo folder:  cd MsaPhylo && python ../pipeline/run_attention_trees.py X N
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, os.getcwd())  # MsaPhylo's code/ package, not the standard-library module of the same name
import code.extracting as extracting
import code.build_plm_tree as build_plm_tree
from code.extracting import Extractor
from code.build_plm_tree import PlmTree

X, N = sys.argv[1], sys.argv[2]
INPUT_DIR = Path(f"data/{X}45_msas_{N}")
ATTN_DIR = Path(f"Attentions/{X}45_{N}")
OUTPUT_DIR = Path(f"./{X}45_attention_trees_{N}")
ATTN_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Redirect the paths in code/params.py so data/Pfam, Attentions/Pfam and Trees/ are left untouched.
extracting.MSA_PATH = build_plm_tree.MSA_PATH = f"{INPUT_DIR}/"
extracting.ATTN_PATH = build_plm_tree.ATTN_PATH = f"{ATTN_DIR}/"

msa_type = "default"
msa_files = sorted(INPUT_DIR.glob("*.fasta"))

print("=" * 80)
print(f"Input directory:  {INPUT_DIR.resolve()}")
print(f"Output directory: {OUTPUT_DIR.resolve()}")
print(f"Found {len(msa_files)} MSA files")
print("=" * 80)

failed = []
# Run MsaPhylo attention-tree inference for all {X}45_{N} files.
for i, msa_file in enumerate(msa_files, start=1):
    # Example: <X>45_supermatrix_<N>.fasta -> <X>45_supermatrix_<N>
    prot_family = msa_file.stem
    tree_dir = OUTPUT_DIR / prot_family
    tree_dir.mkdir(parents=True, exist_ok=True)
    build_plm_tree.TREE_PATH = f"{tree_dir}/"

    print("\n" + "=" * 80)
    print(f"[{i}/{len(msa_files)}] Running {prot_family}")
    print(f"MSA file: {msa_file}")
    print("=" * 80)

    try:
        # Step 1: Extract and save column attentions.
        ex = Extractor(prot_family, msa_type)
        ex.get_col_attention()

        # Step 2: Build NJ trees from every attention head (12 layers x 12 heads).
        plmtree = PlmTree(prot_family, msa_type)
        plmtree.build_attention_tree()

        print(f"Finished {prot_family}")

    except Exception as e:
        print(f"Failed {prot_family}")
        print(f"Error: {e}")
        failed.append(prot_family)
        continue

print(f"\nAll finished. Failed: {len(failed)} {failed}")
