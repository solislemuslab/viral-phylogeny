"""
Circular trees with coloured clades (item b.i).

    python plot_trees.py --tree-dir /path/to/trees --groups /path/to/groups.tsv

With no arguments it runs on the bundled PF00066 test trees. See tree_utils.py
for all shared flags.

Writes <out>/circular_<method>.png (one labelled tree per method) and
<out>/circular_overview.png (all trees side by side).

When a groups file is found (--groups, or see tree_utils.find_groups_file), a
branch takes the group colour when every tip below it belongs to that one
group, otherwise it is grey, and each maximal pure clade gets an arc outside
the tips; group is also encoded by tip marker shape, so identity never rests on
colour alone.  Without one, trees are drawn in plain grey.
draw_circular is reused by the other scripts for their circular panels.
"""
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge

from tree_utils import INK, INK_2, MUTED, SURFACE, load_all, load_groups, parse_args, rooted_all

R0 = 0.03  # radius of the root


def _arc(r, a0, a1, step=math.radians(1.5)):
    k = max(2, int(abs(a1 - a0) / step) + 1)
    t = np.linspace(a0, a1, k)
    return np.column_stack([r * np.cos(t), r * np.sin(t)])


def pure_edge_color(groups):
    """edge_color_fn: group colour for a pure clade, grey otherwise; plain grey
    for every edge when there is no groups file at all."""
    if groups is None:
        return lambda node: MUTED

    def fn(node):
        g = groups.group_of_mask(node.mask)
        return groups.color(g) if g else MUTED
    return fn


def maximal_pure_clades(root, groups):
    """[(node, group)] - the largest clades whose tips are all in one group."""
    out = []

    def rec(nd):
        g = groups.group_of_mask(nd.mask)
        if g and nd.leaf is None:
            out.append((nd, g))
            return
        for c in nd.children:
            rec(c)
    rec(root)
    return out


def draw_circular(ax, root, taxa, *, groups=None, edge_color=None, labels=True,
                  label_size=5.0, label_fn=None, tip_ring=None, gap_deg=4.0, lw=1.0):
    """Draw a ladderized rooted tree (RNode) as a circular cladogram.

    groups     -> tip markers by group and arcs over maximal pure clades
    edge_color -> fn(node) -> colour of the edge above `node`
    tip_ring   -> {taxon index: colour}; a heat-map band outside the tips
    """
    leaves = root.leaves()
    n = len(leaves)
    gap = math.radians(gap_deg)
    span = 2 * math.pi - gap
    ang = {}
    for i, lf in enumerate(leaves):
        ang[lf] = gap / 2 + span * (i + 0.5) / n
    H = root.height

    def radius(nd):
        return R0 + (1 - R0) * (1 - nd.height / H)

    segs, cols = [], []

    def place(nd):
        if nd.children:
            for c in nd.children:
                place(c)
            ang[nd] = (ang[nd.children[0]] + ang[nd.children[-1]]) / 2
    place(root)

    for nd in root.walk():
        if not nd.children:
            continue
        r = radius(nd)
        if nd is not root:
            col = edge_color(nd) if edge_color else MUTED
            a0, a1 = ang[nd.children[0]], ang[nd.children[-1]]
            segs.append(_arc(r, a0, a1))
            cols.append(col)
        for c in nd.children:
            rc = radius(c)
            a = ang[c]
            segs.append(np.array([[r * math.cos(a), r * math.sin(a)],
                                  [rc * math.cos(a), rc * math.sin(a)]]))
            cols.append(edge_color(c) if edge_color else MUTED)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lw, capstyle="round"))

    # tip markers (only when a groups file was found)
    if groups is not None:
        for g in groups.order:
            pts = [(math.cos(ang[l]) * 1.012, math.sin(ang[l]) * 1.012)
                   for l in leaves if groups.broad[taxa[l.leaf]] == g]
            if pts:
                xs, ys = zip(*pts)
                ax.scatter(xs, ys, s=11, marker=groups.marker(g), color=groups.color(g),
                           edgecolors=SURFACE, linewidths=0.4, zorder=3)
        for nd, g in maximal_pure_clades(root, groups):
            ls = nd.leaves()
            a0 = min(ang[l] for l in ls) - span / n * 0.36
            a1 = max(ang[l] for l in ls) + span / n * 0.36
            xy = _arc(1.045, a0, a1)
            ax.plot(xy[:, 0], xy[:, 1], color=groups.color(g), lw=2.6, solid_capstyle="butt")

    if tip_ring is not None:
        half = math.degrees(span / n) / 2
        for lf in leaves:
            if lf.leaf in tip_ring:
                d = math.degrees(ang[lf])
                ax.add_patch(Wedge((0, 0), 1.10, d - half, d + half, width=0.05,
                                   facecolor=tip_ring[lf.leaf], edgecolor=SURFACE, lw=0.4))

    if labels:
        r = 1.13 if tip_ring is not None else 1.075
        for lf in leaves:
            a = ang[lf]
            deg = math.degrees(a)
            flip = math.cos(a) < 0
            ax.text(r * math.cos(a), r * math.sin(a),
                    label_fn(lf.leaf) if label_fn else taxa[lf.leaf],
                    rotation=deg + 180 if flip else deg, rotation_mode="anchor",
                    ha="right" if flip else "left", va="center",
                    fontsize=label_size, color=INK)
    lim = 1.62 if labels else 1.16
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    return ang


def group_legend(ax, groups, **kw):
    """Legend of the groups, or None when there is no groups file."""
    if groups is None:
        return None
    handles = [Line2D([], [], marker=groups.marker(g), color=groups.color(g), linestyle="",
                      markersize=7, label=g) for g in groups.order]
    handles.append(Line2D([], [], color=MUTED, lw=1.5, label="mixed clade"))
    return ax.legend(handles=handles, frameon=False, labelcolor=INK_2, fontsize=9, **kw)


def plot_single(name, root, tree, groups, reference, title, path):
    fig, ax = plt.subplots(figsize=(11, 11), facecolor=SURFACE)
    draw_circular(ax, root, tree.taxa, groups=groups, edge_color=pure_edge_color(groups),
                  label_size=5.2)
    group_legend(ax, groups, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    heading = f"{name}: {title} tree" + (", clades coloured by taxonomic group" if groups else "")
    fig.text(0.02, 0.985, heading, ha="left", va="top", fontsize=12, color=INK)
    fig.text(0.02, 0.962,
             f"{len(tree.splits)} internal splits ({tree.resolution:.0%} resolved); "
             f"rooted at the split closest to the {reference} midpoint root; branch lengths ignored"
             + ("" if groups else " (no taxon-group file found; pass --groups to colour clades)"),
             ha="left", va="top", fontsize=8.5, color=INK_2)
    fig.savefig(path, dpi=170, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def plot_overview(rooted, trees, groups, reference, title, path):
    n = len(rooted)
    ncols = 3 if n > 4 else n
    nrows = math.ceil((n + 1) / ncols)  # +1 panel reserved for the legend
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.3 * ncols, 3.6 * nrows), facecolor=SURFACE)
    axes = list(np.atleast_1d(axes).flat)  # a list, not zip()-over-an-iterator: see below
    for ax, (name, root) in zip(axes[:n], rooted.items()):
        draw_circular(ax, root, trees[name].taxa, groups=groups,
                      edge_color=pure_edge_color(groups), labels=False, lw=0.9)
        ax.set_title(f"{name}  ({len(trees[name].splits)} splits)", loc="left",
                     fontsize=11, color=INK)
    lax = axes[n]
    lax.axis("off")
    if groups is not None:
        group_legend(lax, groups, loc="center left")
        caption = ("Branch and arc colour = taxonomic group\n"
                  "(grey = clade mixes groups).\nAll trees rooted on the split\n"
                  f"closest to the {reference} midpoint root;\nbranch lengths ignored.")
    else:
        caption = ("No taxon-group file found (pass --groups\nto colour clades by group).\n"
                  "All trees rooted on the split closest\n"
                  f"to the {reference} midpoint root;\nbranch lengths ignored.")
    lax.text(0.0, 0.18, caption, transform=lax.transAxes, fontsize=9, color=INK_2, va="top")
    for ax in axes[n + 1:]:  # any remaining unused panels
        ax.axis("off")
    subtitle = "how well do clades follow taxonomy?" if groups is not None else "topology only"
    fig.suptitle(f"{title} trees from {n} methods: {subtitle}",
                 x=0.02, y=0.995, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    groups = load_groups(trees.taxa, args.groups, args.tree_dir, args.title)
    rooted = rooted_all(trees)
    for name, root in rooted.items():
        safe = name.replace(" ", "_")
        plot_single(name, root, trees[name], groups, trees.reference, args.title,
                   args.out / f"circular_{safe}.png")
    plot_overview(rooted, trees, groups, trees.reference, args.title,
                 args.out / "circular_overview.png")
    print("wrote", len(rooted), "circular_<method>.png and circular_overview.png in", args.out)


if __name__ == "__main__":
    main()
