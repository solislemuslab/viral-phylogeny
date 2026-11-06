"""
Co-phylogeny / tanglegrams (item b.ii).

    python plot_cophylo.py --tree-dir /path/to/trees --groups /path/to/groups.tsv

With no arguments it runs on the bundled PF00066 test trees. See tree_utils.py
for all shared flags.

Writes <out>/cophylo_<reference>_vs_<method>.png for the reference tree against
each other method, and <out>/cophylo_crossings.csv with the crossing count for
every pair of methods.

The reference tree keeps its ladderized order.  The other tree has its children
re-ordered (topology unchanged) so each node sits near the mean position of its
tips in the reference, which removes avoidable crossings.  Lines join the same
taxon in both trees, coloured by taxonomic group when a groups file was found
(plain otherwise); the crossing count measures how differently the two methods
arrange the taxa.
"""
import csv
import itertools

import matplotlib.pyplot as plt

from plot_trees import group_legend, pure_edge_color
from tree_utils import (INK, INK_2, SURFACE, load_all, load_groups, parse_args,
                        rf_distance, rooted_all)

TREE_W = 3.8      # inches per tree
LABEL_FS = 5.0


def leaf_positions(root):
    return {lf.leaf: i for i, lf in enumerate(root.leaves())}


def untangle(root, ref_pos):
    """Reorder children so each clade sits at the mean reference position of its tips."""
    def key(nd):
        ls = nd.leaves()
        return sum(ref_pos[l.leaf] for l in ls) / len(ls)

    def rec(nd):
        for c in nd.children:
            rec(c)
        nd.children.sort(key=key)
    rec(root)


def crossings(pos_a, pos_b):
    """Number of taxon pairs whose order differs between the two leaf orders."""
    taxa = list(pos_a)
    return sum((pos_a[i] - pos_a[j]) * (pos_b[i] - pos_b[j]) < 0
               for i, j in itertools.combinations(taxa, 2))


def draw_rect(ax, root, ypos, x_root, x_tip, groups):
    """Rectangular cladogram; tips aligned at x_tip, root at x_root."""
    H = root.height
    edge_color = pure_edge_color(groups)
    ys = {}

    def x_of(nd):
        return x_root + (x_tip - x_root) * (1 - nd.height / H)

    def rec(nd):
        if nd.leaf is not None:
            ys[nd] = ypos[nd.leaf]
            return
        for c in nd.children:
            rec(c)
        ys[nd] = (ys[nd.children[0]] + ys[nd.children[-1]]) / 2
    rec(root)
    for nd in root.walk():
        if not nd.children:
            continue
        x = x_of(nd)
        ax.plot([x, x], [ys[nd.children[0]], ys[nd.children[-1]]],
                color=edge_color(nd), lw=0.9, solid_capstyle="round")
        for c in nd.children:
            ax.plot([x, x_of(c)], [ys[c], ys[c]], color=edge_color(c), lw=0.9,
                    solid_capstyle="round")


def plot_pair(name_a, name_b, root_a, root_b, tree_a, tree_b, groups, path):
    taxa = tree_a.taxa
    n = len(taxa)
    pos_a, pos_b = leaf_positions(root_a), leaf_positions(root_b)
    cross = crossings(pos_a, pos_b)
    rf, mx = rf_distance(tree_a, tree_b)

    lab_w = max(len(t) for t in taxa) * LABEL_FS * 0.58 / 72
    gap = 2 * lab_w + 3.2
    width = 2 * TREE_W + gap + 0.6
    height = 15.5
    fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 0.94])
    ax.set_xlim(0, width)
    ax.set_ylim(-4.5, n + 0.2)
    ax.axis("off")

    xl_tip = 0.3 + TREE_W
    xr_tip = xl_tip + gap
    ya = {i: n - 1 - p for i, p in pos_a.items()}
    yb = {i: n - 1 - p for i, p in pos_b.items()}
    draw_rect(ax, root_a, ya, 0.3, xl_tip, groups)
    draw_rect(ax, root_b, yb, xr_tip + TREE_W, xr_tip, groups)

    for i, t in enumerate(taxa):
        col = groups.color(groups.broad[t]) if groups is not None else INK_2
        ax.text(xl_tip + 0.14, ya[i], t, ha="left", va="center", fontsize=LABEL_FS, color=INK)
        ax.text(xr_tip - 0.14, yb[i], t, ha="right", va="center", fontsize=LABEL_FS, color=INK)
        for x, y, off in ((xl_tip, ya[i], 0.05), (xr_tip, yb[i], -0.05)):
            marker = groups.marker(groups.broad[t]) if groups is not None else "o"
            ax.scatter([x + off], [y], s=9, marker=marker, color=col,
                       edgecolors=SURFACE, linewidths=0.3, zorder=3)
        ax.plot([xl_tip + 0.14 + lab_w + 0.05, xr_tip - 0.14 - lab_w - 0.05],
                [ya[i], yb[i]], color=col, lw=0.8, alpha=0.85, zorder=1)

    if groups is not None:
        group_legend(ax, groups, loc="lower center", ncol=len(groups.order) + 1,
                     bbox_to_anchor=(0.5, 0.0))
    fig.text(0.015, 0.985, f"{name_a}  vs  {name_b}", ha="left", va="top", fontsize=13,
             color=INK)
    fig.text(0.015, 0.962,
             f"{cross} crossing taxon pairs of {n * (n - 1) // 2}; "
             f"normalized RF {rf / mx:.2f}. Lines join the same taxon"
             + (", coloured by taxonomic group; " if groups is not None else "; ")
             + f"{name_b} re-ordered to match {name_a}.",
             ha="left", va="top", fontsize=8.5, color=INK_2)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return cross


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    groups = load_groups(trees.taxa, args.groups, args.tree_dir, args.title)
    names = list(trees)
    reference = trees.reference

    # reference vs each other method (reference order stays fixed across figures)
    ref_root = rooted_all(trees)[reference]
    ref_pos = leaf_positions(ref_root)
    for name in names:
        if name == reference:
            continue
        other = rooted_all(trees)[name]
        untangle(other, ref_pos)
        path = args.out / f"cophylo_{reference}_vs_{name}.png".replace(" ", "_")
        c = plot_pair(reference, name, ref_root, other, trees[reference], trees[name],
                      groups, path)
        print(f"{reference} vs {name}: {c} crossings -> {path.name}")

    # crossing counts for every pair (second tree re-ordered to match the first)
    rows = []
    for a, b in itertools.combinations(names, 2):
        ra, rb = rooted_all(trees)[a], rooted_all(trees)[b]
        untangle(rb, leaf_positions(ra))
        rows.append((a, b, crossings(leaf_positions(ra), leaf_positions(rb))))
    with open(args.out / "cophylo_crossings.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tree_a", "tree_b", "crossings"])
        w.writerows(rows)
    print("wrote cophylo_crossings.csv")


if __name__ == "__main__":
    main()
