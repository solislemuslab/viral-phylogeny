"""
Write Comparison/compare_simple_<X>_<label>.R for one group of result trees.

    python pipeline/make_compare_R.py X label METHODS
    e.g. python pipeline/make_compare_R.py AS 113 IQ,PF,CM

label is the group name used in <X>45_<label>_result_trees (a number of random loci, or e.g. short_26).
METHODS is a comma list from IQ (IQ-TREE), PF (Phyloformer), CM (CMAPLE), MA (MSA-T attention), ME (MSA-T embedding).
The script is pipeline/compare_template.R (same plots as Comparison/compare_simple_OB_603.R) with the dataset's
outgroup (pipeline/config.sh), result-tree files, reference trees and clades (pipeline/clades/<X>.R) filled in.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIPE = ROOT / "pipeline"
X, label, methods = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
outgroup = re.search(rf'^\s*{X}\)\s*OUTGROUP=(\w+)', (PIPE / "config.sh").read_text(), re.M).group(1)

FILES = {"IQ": ("IQ-TREE", "{X}45_supermatrix_{n}_IQTree.treefile"),
         "PF": ("Phyloformer", "{X}45_supermatrix_{n}_Phyloformer.nwk"),
         "CM": ("CMAPLE", "{X}45_supermatrix_{n}_CMAPLE.treefile"),
         "MA": ("MSA-T attention", "{X}45_supermatrix_{n}_0_4_MSATransformer_AttentionTree.nwk"),
         "ME": ("MSA-T embedding", "{X}45_{n}_2_embedding_tree.nwk")}
lines = [f'           {chr(34) + FILES[m][0] + chr(34):<17} = "{FILES[m][1].format(X=X, n=label)}"' for m in methods]
lines[0] = "files <- c(" + lines[0].lstrip()
header = (f'outgroup <- "{outgroup}" # root of {X}-astral.tre; used for every plot (RF/quartets stay unrooted)\n'
          f'tree_dir <- "../{X}45_{label}_result_trees"\n' + ",\n".join(lines) + ")\n"
          f'ref_dir <- "../{X}45_reference_trees"\n'
          f'ref_files <- c("{X}-astral"      = "{X}-astral.tre",\n'
          f'               "{X}-ncbi-astral" = "{X}-ncbi-astral.tre",\n'
          f'               "{X}-ncbi"        = "{X}-ncbi.tre")')
s = (PIPE / "compare_template.R").read_text().split("\n", 1)[1]
s = s.replace("@@HEADER@@", header).replace("@@CLADES@@", (PIPE / f"clades/{X}.R").read_text())
s = s.replace("@@X@@", X).replace("@@LABEL@@", label)
out = ROOT / f"Comparison/compare_simple_{X}_{label}.R"
out.write_text(s)
print("wrote", out)
