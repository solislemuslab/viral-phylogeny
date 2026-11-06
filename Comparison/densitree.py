"""
DensiTree on a reference tree (item b.iv).

    python densitree.py --tree-dir /path/to/trees --groups /path/to/groups.tsv

With no arguments it runs on the bundled PF00066 test trees. See tree_utils.py
for all shared flags.

Writes <out>/densitree.png (every tree overlaid, one colour per method) and
<out>/densitree_by_method.png (one panel per method, each over the reference in
grey, so no colour has to be told apart from several others).

Tip order is the reference tree's, so where the methods agree their lines run
together; where they disagree the lines fan out.  Every tree is rooted with the
same convention as the circular trees and drawn as a cladogram scaled to its own
root (branch lengths ignored).
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

from tree_utils import INK, INK_2, MUTED, NEUTRAL, SURFACE, load_all, load_groups, parse_args, rooted_all

N_PTS = 28
T = np.linspace(0, 1, N_PTS)
SMOOTH = T * T * (3 - 2 * T)  # curve leaves the parent and enters the child horizontally


def tree_segments(root, ypos):
    """Polylines for every parent->child edge; x = depth (root 0, tips 1)."""
    H = root.height
    x, y = {}, {}

    def rec(nd):
        x[nd] = 1 - nd.height / H
        if nd.leaf is not None:
            y[nd] = ypos[nd.leaf]
            return
        for c in nd.children:
            rec(c)
        y[nd] = sum(y[c] for c in nd.children) / len(nd.children)
    rec(root)
    segs = []
    for nd in root.walk():
        for c in nd.children:
            xs = x[nd] + (x[c] - x[nd]) * T
            ys = y[nd] + (y[c] - y[nd]) * SMOOTH
            segs.append(np.column_stack([xs, ys]))
    return segs


def draw_tip_labels(ax, taxa, ypos, groups, size):
    for i, t in enumerate(taxa):
        color = groups.color(groups.broad[t]) if groups is not None else MUTED
        marker = groups.marker(groups.broad[t]) if groups is not None else "o"
        ax.scatter([1.012], [ypos[i]], s=size * 2.2, marker=marker, color=color,
                   edgecolors=SURFACE, linewidths=0.3, zorder=3)
        ax.text(1.028, ypos[i], t, ha="left", va="center", fontsize=size, color=INK)


def plot_overlay(rooted, taxa, ypos, groups, colors, reference, title, path):
    n = len(taxa)
    fig, ax = plt.subplots(figsize=(11.5, 17), facecolor=SURFACE)
    ax.add_collection(LineCollection(tree_segments(rooted[reference], ypos), colors=NEUTRAL,
                                     linewidths=3.4, capstyle="round", zorder=1))
    for name, root in rooted.items():
        ax.add_collection(LineCollection(tree_segments(root, ypos),
                                         colors=colors[name], linewidths=0.9,
                                         alpha=0.55, capstyle="round", zorder=2))
    draw_tip_labels(ax, taxa, ypos, groups, 5.0)
    ax.set_xlim(-0.02, 1.36)
    ax.set_ylim(n + 1, -3.5)
    ax.axis("off")
    handles = [Line2D([], [], color=colors[m], lw=2.2, label=m) for m in rooted]
    handles.append(Line2D([], [], color=NEUTRAL, lw=4, label=f"{reference} topology (backdrop)"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False,
              labelcolor=INK_2, fontsize=9)
    fig.text(0.02, 0.99, f"DensiTree of the {len(rooted)} {title} trees on the {reference} tip order",
             ha="left", va="top", fontsize=12, color=INK)
    fig.text(0.02, 0.972,
             "root on the left, tips on the right; every tree drawn as a cladogram scaled to its own root",
             ha="left", va="top", fontsize=8.5, color=INK_2)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def plot_panels(rooted, taxa, ypos, groups, colors, reference, path):
    n = len(taxa)
    fig, axes = plt.subplots(1, len(rooted), figsize=(4.6 * len(rooted), 9), facecolor=SURFACE)
    axes = np.atleast_1d(axes)
    ref_segs = tree_segments(rooted[reference], ypos)
    for ax, (name, root) in zip(axes, rooted.items()):
        ax.add_collection(LineCollection(ref_segs, colors=NEUTRAL, linewidths=2.6,
                                         capstyle="round", zorder=1))
        if name != reference:
            ax.add_collection(LineCollection(tree_segments(root, ypos),
                                             colors=colors[name], linewidths=0.9,
                                             alpha=0.8, capstyle="round", zorder=2))
        else:
            ax.add_collection(LineCollection(ref_segs, colors=colors[name],
                                             linewidths=0.9, alpha=0.9, zorder=2))
        for i, t in enumerate(taxa):
            color = groups.color(groups.broad[t]) if groups is not None else MUTED
            marker = groups.marker(groups.broad[t]) if groups is not None else "o"
            ax.scatter([1.02], [ypos[i]], s=5, marker=marker, color=color, linewidths=0, zorder=3)
        ax.set_xlim(-0.02, 1.06)
        ax.set_ylim(n, -1)
        ax.axis("off")
        ax.set_title(name, loc="left", fontsize=11, color=INK)
    fig.suptitle(f"Each method against the {reference} topology (grey); tips in {reference} order",
                 x=0.01, y=0.995, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    taxa = trees.taxa
    groups = load_groups(taxa, args.groups, args.tree_dir, args.title)
    rooted = rooted_all(trees)
    order = rooted[trees.reference].leaves()
    ypos = {lf.leaf: i for i, lf in enumerate(order)}
    plot_overlay(rooted, taxa, ypos, groups, trees.colors, trees.reference, args.title,
                args.out / "densitree.png")
    plot_panels(rooted, taxa, ypos, groups, trees.colors, trees.reference,
               args.out / "densitree_by_method.png")
    print("wrote densitree.png and densitree_by_method.png in", args.out)


if __name__ == "__main__":
    main()
