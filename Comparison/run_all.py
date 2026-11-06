"""
Run every comparison in order, for one tree set, and write all figures/tables to
one output directory.

    conda activate phylo        # needs dendropy, numpy, matplotlib
    python run_all.py --tree-dir /path/to/trees --groups /path/to/groups.tsv

With no arguments it runs on the bundled PF00066 test trees
(../PF00066_result _trees, groups from PF00066_groups.tsv; see
PF00066_result _trees/trees.tsv for how its five trees are named).  See
tree_utils.py for --tree-dir/--reference/--groups/--out/--title, which apply to
every step below.
"""
import compare
import consensus
import densitree
import monophyly
import plot_cophylo
import plot_trees
import quartets
import reference_clades
import tree_utils

STEPS = [
    ("RF matrix (iii)", compare.main),
    ("quartet consistency (iii)", quartets.main),
    ("circular trees (i)", plot_trees.main),
    ("co-phylogeny (ii)", plot_cophylo.main),
    ("densitree (iv)", densitree.main),
    ("reference clades (v)", reference_clades.main),
    ("consensus (vi)", consensus.main),
    ("monophyly (vii)", monophyly.main),
]

if __name__ == "__main__":
    args = tree_utils.parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"dataset: {args.title}  ({args.tree_dir})")
    print(f"writing to {args.out}")
    for title, step in STEPS:
        print(f"\n=== {title}")
        step(args)
