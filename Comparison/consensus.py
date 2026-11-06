"""
Consensus of every tree in a tree set (item b.vi).

    python consensus.py --tree-dir /path/to/trees --groups /path/to/groups.tsv

With no arguments it runs on the bundled PF00066 test trees. See tree_utils.py
for all shared flags.

Split frequencies are counted over the unrooted trees, then three consensus
trees are built from them:
  strict    splits found in every tree
  majority  splits found in more than half the trees
  greedy    splits added in order of frequency whenever compatible with those already
            kept (extended majority rule) - the most resolved of the three
A consensus is unresolved (polytomies) wherever the methods disagree.

Writes to <out>/:
  consensus_strict.nwk / consensus_majority.nwk / consensus_greedy.nwk
      internal node labels are the number of trees containing that clade;
      rooted on the split closest to the reference's midpoint root, as in the drawings
  consensus_splits.csv       every split seen in >= 2 trees, with frequency and membership
  consensus_overview.png     the three consensus trees side by side
  consensus_majority.png     majority-rule tree with tip labels
"""
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from plot_trees import draw_circular
from tree_utils import (INK, INK_2, NEUTRAL, SURFACE, best_match_root, blues, compatible,
                        ladderize, load_all, load_groups, parse_args, popcount, root_masks,
                        rooted_tree, to_newick)


def split_frequencies(trees):
    """{split: [methods that contain it]}"""
    freq = defaultdict(list)
    for name, t in trees.items():
        for s in t.splits:
            freq[s].append(name)
    return freq


def build_consensus(freq, n_trees):
    strict = {s for s, m in freq.items() if len(m) == n_trees}
    majority = {s for s, m in freq.items() if len(m) * 2 > n_trees}
    greedy = []
    for s in sorted(freq, key=lambda s: (-len(freq[s]), -popcount(s), s)):
        if all(compatible(s, t) for t in greedy):
            greedy.append(s)
    return {"strict": strict, "majority": majority, "greedy": set(greedy)}


def rooted_consensus(splits, n, ref_root):
    root_mask = best_match_root(splits, ref_root, n)
    return ladderize(rooted_tree(splits, n, root_mask))


def plot_tree(ax, root, freq, n, taxa, steps, labels, label_size=4.8):
    def edge_color(nd):
        if nd.leaf is not None:
            return NEUTRAL
        return steps[len(freq[nd.split(n)]) - 1]
    draw_circular(ax, root, taxa, edge_color=edge_color, labels=labels,
                  label_size=label_size, lw=1.4 if labels else 1.2)


def freq_legend(ax, n_trees, steps, **kw):
    handles = [Line2D([], [], color=NEUTRAL, lw=2.5, label="tip edge")]
    handles += [Line2D([], [], color=steps[k - 1], lw=2.5, label=f"clade in {k} of {n_trees} trees")
                for k in range(n_trees, 0, -1)]
    return ax.legend(handles=handles, frameon=False, labelcolor=INK_2, fontsize=9, **kw)


def largest_polytomy(root):
    return max((len(nd.children) for nd in root.walk() if nd.children), default=0)


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    reference = trees.reference
    ref = trees[reference]
    n, taxa, n_trees = ref.n, trees.taxa, len(trees)
    groups = load_groups(taxa, args.groups, args.tree_dir, args.title)
    ref_root = root_masks(trees)[reference]
    freq = split_frequencies(trees)
    cons = build_consensus(freq, n_trees)
    steps = blues(n_trees)

    roots = {}
    for kind, splits in cons.items():
        root = rooted_consensus(splits, n, ref_root)
        roots[kind] = root
        (args.out / f"consensus_{kind}.nwk").write_text(
            to_newick(root, taxa, lambda nd: len(freq[nd.split(n)])) + ";\n")
        print(f"{kind:9s} {len(splits):3d} of {n - 3} splits ({len(splits) / (n - 3):.0%} resolved), "
              f"largest polytomy {largest_polytomy(root)} branches")

    with open(args.out / "consensus_splits.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["size", "n_trees", "in_strict", "in_majority", "in_greedy", "composition",
                    "methods", "members"])
        for s in sorted(freq, key=lambda s: (-len(freq[s]), -popcount(s), s)):
            if len(freq[s]) < 2:
                continue
            side = s if popcount(s) <= n - popcount(s) else ((1 << n) - 1) ^ s
            members = [taxa[i] for i in range(n) if (side >> i) & 1]
            if groups is not None:
                comp = defaultdict(int)
                for t in members:
                    comp[groups.broad[t]] += 1
                top = max(comp, key=comp.get)
                composition = f"{top} {comp[top]}/{len(members)}"
            else:
                composition = ""
            w.writerow([len(members), len(freq[s]), int(s in cons["strict"]),
                        int(s in cons["majority"]), int(s in cons["greedy"]),
                        composition, ";".join(freq[s]), ";".join(members)])

    fig, axes = plt.subplots(1, 3, figsize=(17, 6.4), facecolor=SURFACE)
    for ax, (kind, root) in zip(axes, roots.items()):
        plot_tree(ax, root, freq, n, taxa, steps, labels=False)
        ax.set_title(f"{kind}: {len(cons[kind])} of {n - 3} splits", loc="left", fontsize=11,
                     color=INK)
    freq_legend(fig, n_trees, steps, loc="lower center", ncol=n_trees + 1)
    fig.suptitle(f"Consensus of the {n_trees} {args.title} trees; "
                "branches coloured by how many trees contain the clade",
                x=0.01, y=0.995, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(args.out / "consensus_overview.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 11), facecolor=SURFACE)
    plot_tree(ax, roots["majority"], freq, n, taxa, steps, labels=True)
    freq_legend(ax, n_trees, steps, loc="upper left", bbox_to_anchor=(0.0, 0.99))
    fig.text(0.02, 0.985, f"Majority-rule consensus of the {n_trees} {args.title} trees",
             ha="left", va="top", fontsize=12, color=INK)
    fig.text(0.02, 0.962, f"{len(cons['majority'])} of {n - 3} possible splits kept "
             f"(clades present in at least {n_trees // 2 + 1} of {n_trees} trees); "
             "polytomies mark disagreement", ha="left", va="top", fontsize=8.5, color=INK_2)
    fig.savefig(args.out / "consensus_majority.png", dpi=170, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote consensus_*.nwk, consensus_splits.csv, consensus_overview.png, consensus_majority.png")


if __name__ == "__main__":
    main()
