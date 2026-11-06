"""
Which reference clades do the other methods recover? (item b.v)

    python reference_clades.py --tree-dir /path/to/trees --groups /path/to/groups.tsv

With no arguments it runs on the bundled PF00066 test trees. See tree_utils.py
for all shared flags.

Every internal edge of the reference tree defines a clade.  A clade counts as
recovered by another method when that tree has the same split (unrooted, so the
answer does not depend on where either tree is rooted; a polytomy that swallows
the clade counts as not recovered).

Writes to <out>/:
  reference_clades.csv          one row per reference clade: size, tip composition
                                (when a groups file was found), reference support,
                                and present/absent in each method
  reference_clades_summary.csv  per method: reference splits recovered (recall) and
                                the share of the method's own splits found in the reference
  reference_clades.png          reference tree with each branch coloured by how many
                                of the other methods recover it
"""
import csv
from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from plot_trees import draw_circular
from tree_utils import (INK, INK_2, NEUTRAL, SURFACE, blues, load_all, load_groups, parse_args,
                        rooted_all)


def clade_table(trees, groups, ref_root, reference):
    """One row per reference split, described by the clade side away from the root."""
    ref = trees[reference]
    others = [m for m in trees if m != reference]
    n = ref.n
    rows = {}
    for nd in ref_root.walk():
        if nd.leaf is not None or nd is ref_root:
            continue
        s = nd.split(n)
        if s in rows:  # the two root children share one split
            continue
        if groups is not None:
            comp = Counter(groups.broad[ref.taxa[l.leaf]] for l in nd.leaves())
            top, k = comp.most_common(1)[0]
            composition = f"{top} {k}/{nd.size}"
        else:
            composition = ""
        present = {m: s in trees[m].splits for m in others}
        rows[s] = {
            "split": s,
            "size": nd.size,
            "support": ref.splits[s][1] or "",
            "composition": composition,
            "n_methods": sum(present.values()),
            "present": present,
            "members": ";".join(ref.taxa[l.leaf] for l in nd.leaves()),
        }
    return sorted(rows.values(), key=lambda r: (-r["size"], r["split"])), others


def method_summary(trees, reference):
    ref = trees[reference].splits.keys()
    out = []
    for m, t in trees.items():
        if m == reference:
            continue
        shared = len(ref & t.splits.keys())
        out.append((m, shared, len(ref), shared / len(ref), shared / len(t.splits)))
    return out


def plot_recovery(ref_root, ref, rows, summary, colors, reference, path):
    count = {r["split"]: r["n_methods"] for r in rows}
    n = ref.n
    n_oth = len(summary)
    steps = blues(n_oth)

    def edge_color(nd):
        if nd.leaf is not None:
            return NEUTRAL
        k = count[nd.split(n)]
        return NEUTRAL if k == 0 else steps[k - 1]

    fig = plt.figure(figsize=(14, 11), facecolor=SURFACE)
    ax = fig.add_axes([0.0, 0.0, 0.76, 0.94])
    draw_circular(ax, ref_root, ref.taxa, edge_color=edge_color, label_size=4.8, lw=1.5)
    handles = [Line2D([], [], color=NEUTRAL, lw=2.5, label="in no other method")]
    handles += [Line2D([], [], color=steps[k - 1], lw=2.5,
                       label=f"recovered by {k} of {n_oth}") for k in range(1, n_oth + 1)]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 0.99), frameon=False,
              labelcolor=INK_2, fontsize=9)

    bx = fig.add_axes([0.87, 0.60, 0.11, 0.16], facecolor=SURFACE)
    names = [s[0] for s in summary][::-1]
    vals = [s[3] * 100 for s in summary][::-1]
    bx.barh(range(len(names)), vals, color=[colors[m] for m in names], height=0.6)
    for i, (v, s) in enumerate(zip(vals, summary[::-1])):
        bx.text(v + 1, i, f"{v:.0f}%  ({s[1]}/{s[2]})", va="center", fontsize=8, color=INK)
    bx.set_yticks(range(len(names)), names, fontsize=8, color=INK)
    bx.set_xlim(0, 100)
    bx.set_xticks([])
    bx.tick_params(length=0)
    for side in bx.spines.values():
        side.set_visible(False)
    bx.set_title(f"Share of {reference} clades\nfound by each method", loc="left",
                 fontsize=8.5, color=INK)

    fig.text(0.02, 0.985, f"Which {reference} clades do the other methods recover?",
             ha="left", va="top", fontsize=12, color=INK)
    fig.text(0.02, 0.962, "branch = clade above it; colour = number of other methods with exactly "
             "that split (unrooted); tip edges not counted",
             ha="left", va="top", fontsize=8.5, color=INK_2)
    fig.savefig(path, dpi=170, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    reference = trees.reference
    ref = trees[reference]
    groups = load_groups(trees.taxa, args.groups, args.tree_dir, args.title)
    ref_root = rooted_all(trees)[reference]
    rows, others = clade_table(trees, groups, ref_root, reference)

    with open(args.out / "reference_clades.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["clade", "size", "composition", f"{reference}_support", "n_methods"]
                   + others + ["members"])
        for i, r in enumerate(rows, 1):
            w.writerow([i, r["size"], r["composition"], r["support"], r["n_methods"]]
                       + ["yes" if r["present"][m] else "no" for m in others] + [r["members"]])
    summary = method_summary(trees, reference)
    with open(args.out / "reference_clades_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "shared_splits", "reference_splits", "recall", "precision"])
        for m, shared, tot, rec, prec in summary:
            w.writerow([m, shared, tot, f"{rec:.4f}", f"{prec:.4f}"])
            print(f"{m:16s} recovers {shared:3d}/{tot} reference clades ({rec:.0%})")
    dist = Counter(r["n_methods"] for r in rows)
    print("reference clades by number of other methods recovering them:",
          {k: dist[k] for k in sorted(dist)})
    plot_recovery(ref_root, ref, rows, summary, trees.colors, reference,
                 args.out / "reference_clades.png")
    print("wrote reference_clades*.csv and reference_clades.png in", args.out)


if __name__ == "__main__":
    main()
