"""
Quartet consistency across every tree in a tree set (item b.iii, points 1 and 2).

    python quartets.py --tree-dir /path/to/trees

With no arguments it runs on the bundled PF00066 test trees. See tree_utils.py
for all shared flags.

For every 4-taxon subset a tree either resolves it (ab|cd, ac|bd or ad|bc) or
leaves it as a star.  A quartet is a star exactly when no retained internal edge
separates the four taxa, so polytomies (edges of length ~0, see
tree_utils.COLLAPSE_TOL) are counted as "unresolved", never as a disagreement.
Resolved-in-both is the only thing compared between two trees.

Writes to <out>/:
  quartet_agreement.csv / .png   pairwise: of quartets resolved in both trees,
                                 the share with the same topology
  quartet_levels.png             how many trees agree on each quartet
  quartet_taxon_stability.png    per taxon: share of its quartets on which all
                                 trees agree (circular reference tree)
  quartet_taxon_stability.csv    the same numbers as a table
  quartet_summary.csv            resolved-quartet share per tree
  quartets_common.npz            the quartets every tree resolves identically
                                 (taxon indices a<b<c<d, topology 1/2/3 =
                                 ab|cd, ac|bd, ad|bc; `taxa` gives the names)
"""
import csv
import itertools
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from compare import BLUES, plot_heatmap
from plot_trees import draw_circular
from tree_utils import INK, INK_2, SURFACE, load_all, parse_args, rooted_all


def all_quartets(n):
    """(Q, 4) int16 array of every a<b<c<d."""
    q = math.comb(n, 4)
    flat = np.fromiter(itertools.chain.from_iterable(itertools.combinations(range(n), 4)),
                       dtype=np.int16, count=4 * q)
    return flat.reshape(q, 4)


def path_lengths(tree):
    """Leaf-to-leaf path length in retained edges (collapsed edges count 0)."""
    n = tree.n
    d = np.full((n, n), 2, dtype=np.int16)  # two pendant edges
    np.fill_diagonal(d, 0)
    for s in tree.splits:
        side = np.array([(s >> i) & 1 for i in range(n)], dtype=bool)
        d += (side[:, None] != side[None, :]).astype(np.int16)
    return d


def quartet_codes(d, quartets):
    """0 = star, 1 = ab|cd, 2 = ac|bd, 3 = ad|bc (four-point condition)."""
    a, b, c, e = quartets.T
    s1 = d[a, b] + d[c, e]
    s2 = d[a, c] + d[b, e]
    s3 = d[a, e] + d[b, c]
    mn = np.minimum(np.minimum(s1, s2), s3)
    m1, m2, m3 = s1 == mn, s2 == mn, s3 == mn
    unique = (m1.astype(np.int8) + m2 + m3) == 1
    code = m1 * 1 + m2 * 2 + m3 * 3
    return np.where(unique, code, 0).astype(np.uint8)


def pairwise_agreement(codes, names):
    k = len(names)
    mat = [[1.0] * k for _ in range(k)]
    both_counts = {}
    for i, j in itertools.combinations(range(k), 2):
        both = (codes[i] > 0) & (codes[j] > 0)
        agree = (codes[i] == codes[j]) & both
        mat[i][j] = mat[j][i] = agree.sum() / both.sum()
        both_counts[(names[i], names[j])] = int(both.sum())
    return mat, both_counts


def plot_levels(counts, total, n_trees, path):
    labels = [f"all {n_trees} agree"] + [f"{k} of {n_trees} agree" for k in range(n_trees - 1, 1, -1)]
    labels.append("no two agree")
    vals = [counts[n_trees - i] if n_trees - i >= 2 else counts[0] + counts[1] for i in range(n_trees)]
    fig, ax = plt.subplots(figsize=(8, 3.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    y = np.arange(len(vals))[::-1]
    ax.barh(y, [v / total * 100 for v in vals], color="#2a78d6", height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v / total * 100 + 0.8, yi, f"{v / total:.1%}  ({v:,})", va="center", fontsize=9,
                color=INK)
    ax.set_yticks(y, labels, fontsize=9.5, color=INK)
    ax.set_xlim(0, max(vals) / total * 100 * 1.35)
    ax.set_xticks([])
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title(f"Agreement among the {n_trees} trees on each of the {total:,} quartets\n"
                 "(a tree that leaves a quartet unresolved does not vote)",
                 loc="left", fontsize=10, color=INK)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def plot_taxon_stability(root, taxa, frac, reference, path):
    vmin, vmax = min(frac), max(frac)
    norm = Normalize(vmin, vmax)
    ring = {i: BLUES(norm(f)) for i, f in enumerate(frac)}
    fig, ax = plt.subplots(figsize=(11, 11), facecolor=SURFACE)
    draw_circular(ax, root, taxa, tip_ring=ring, label_size=5.0)
    cax = fig.add_axes([0.08, 0.86, 0.22, 0.014])
    cb = fig.colorbar(ScalarMappable(norm, BLUES), cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.set_label("share of the taxon's quartets on which all trees agree", fontsize=8,
                 color=INK_2)
    cb.ax.tick_params(labelsize=8, length=0, colors=INK_2)
    fig.text(0.02, 0.985, "Where the trees agree on quartets", ha="left", va="top",
             fontsize=12, color=INK)
    fig.text(0.02, 0.962, f"{reference} topology (midpoint-rooted); the outer band is darker "
             "for taxa whose quartets are placed identically by all methods",
             ha="left", va="top", fontsize=8.5, color=INK_2)
    fig.savefig(path, dpi=170, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    names = list(trees)
    reference = trees.reference
    taxa = trees.taxa
    n = len(taxa)
    quartets = all_quartets(n)
    total = len(quartets)
    codes = np.stack([quartet_codes(path_lengths(trees[m]), quartets) for m in names])
    print(f"{total:,} quartets over {n} taxa")

    # per tree: how many quartets are resolved
    with open(args.out / "quartet_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tree", "resolved_quartets", "star_quartets", "resolved_share"])
        for m, c in zip(names, codes):
            res = int((c > 0).sum())
            w.writerow([m, res, total - res, f"{res / total:.5f}"])
            print(f"  {m:16s} resolved {res / total:.2%}")

    # pairwise agreement among mutually resolved quartets
    mat, both = pairwise_agreement(codes, names)
    with open(args.out / "quartet_agreement.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + names)
        for m, row in zip(names, mat):
            w.writerow([m] + [f"{v:.4f}" for v in row])
    agree_title = (f"Quartet agreement, {args.title} trees\n"
                  "share of quartets resolved in both trees that have the same topology")
    plot_heatmap(names, mat, args.out / "quartet_agreement.png", agree_title,
                metric="quartet agreement")
    print("pairwise quartet agreement:")
    for m, row in zip(names, mat):
        print(f"  {m:16s}", *(f"{v:6.3f}" for v in row))

    # how many trees agree on each quartet
    votes = np.stack([(codes == k).sum(axis=0) for k in (1, 2, 3)]).astype(np.uint8)
    top = votes.max(axis=0)
    counts = np.bincount(top, minlength=len(names) + 1)
    plot_levels(counts, total, len(names), args.out / "quartet_levels.png")
    unanimous = top == len(names)
    print(f"quartets all {len(names)} trees resolve identically: "
          f"{unanimous.sum():,} ({unanimous.mean():.2%})")

    # the shared quartets themselves
    idx = np.flatnonzero(unanimous)
    np.savez_compressed(args.out / "quartets_common.npz", quartets=quartets[idx],
                        topology=codes[0][idx], taxa=np.array(taxa))

    # per-taxon stability
    per_taxon_total = math.comb(n - 1, 3)
    hits = np.bincount(quartets[idx].ravel(), minlength=n)
    frac = hits / per_taxon_total
    order = np.argsort(frac)
    with open(args.out / "quartet_taxon_stability.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["taxon", "unanimous_quartets", "share"])
        for i in order[::-1]:
            w.writerow([taxa[i], int(hits[i]), f"{frac[i]:.4f}"])
    print("least stable taxa:", ", ".join(f"{taxa[i]} ({frac[i]:.1%})" for i in order[:5]))
    print("most stable taxa: ", ", ".join(f"{taxa[i]} ({frac[i]:.1%})" for i in order[::-1][:5]))
    root = rooted_all(trees)[reference]
    plot_taxon_stability(root, taxa, frac, reference, args.out / "quartet_taxon_stability.png")
    print("wrote quartet_*.csv / .png and quartets_common.npz in", args.out)


if __name__ == "__main__":
    main()
