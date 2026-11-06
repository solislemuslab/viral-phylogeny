"""
Robinson-Foulds comparison of every tree in a tree set (item b.iii).

    python compare.py --tree-dir /path/to/trees --groups /path/to/groups.tsv

With no arguments it runs on the bundled PF00066 test trees. See tree_utils.py
for all shared flags (--tree-dir/--reference/--groups/--out/--title).

Writes <out>/rf_matrix.csv (normalized RF) and <out>/rf_heatmap.png.
RF is computed on unrooted splits; an edge of length ~0 counts as a polytomy
(tree_utils.COLLAPSE_TOL) and RF is normalized by |splits(a)| + |splits(b)|, so
an unresolved tree is not penalised.

The other comparison scripts (each takes the same flags):
  plot_trees.py        circular trees, coloured clades           (b.i)
  plot_cophylo.py      tanglegrams                                (b.ii)
  quartets.py          quartet consistency, shared quartets       (b.iii)
  densitree.py         DensiTree overlay on the reference         (b.iv)
  reference_clades.py  reference clades present in each tree      (b.v)
  consensus.py         strict / majority / greedy consensus       (b.vi)
  monophyly.py         are groups (taxonomy / region) monophyletic (b.vii)
run_all.py runs all of them on one tree set.
"""
import csv
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from tree_utils import INK, SURFACE, load_all, parse_args, rf_distance

# sequential blue ramp (light -> dark) from the dataviz reference palette
BLUES = LinearSegmentedColormap.from_list("blues", ["#cde2fb", "#3987e5", "#0d366b"])


def compute_rf_matrix(trees):
    """Normalized RF (unrooted) for every pair of trees."""
    names = list(trees)
    mat = [[0.0] * len(names) for _ in names]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i < j:
                rf, max_rf = rf_distance(trees[a], trees[b])
                mat[i][j] = mat[j][i] = rf / max_rf if max_rf else 0.0
    return names, mat


def print_matrix(names, mat):
    width = max(len(n) for n in names)
    print(" " * width, *(f"{i + 1:>6}" for i in range(len(names))))
    for i, (name, row) in enumerate(zip(names, mat)):
        print(f"{name:<{width}}", *(f"{v:6.3f}" for v in row), f"  [{i + 1}]")


def save_csv(names, mat, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + names)
        for name, row in zip(names, mat):
            w.writerow([name] + [f"{v:.4f}" for v in row])


def plot_heatmap(names, mat, path, title, metric="normalized RF"):
    n = len(names)
    off = [mat[i][j] for i in range(n) for j in range(n) if i != j]
    # the off-diagonal values sit in a narrow band, so stretch the scale over it
    vmin = math.floor(min(off) * 20) / 20
    vmax = math.ceil(max(off) * 20) / 20
    shown = np.array(mat, dtype=float)
    np.fill_diagonal(shown, np.nan)  # a tree vs itself carries no information

    fig, ax = plt.subplots(figsize=(8.2, 6.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    cmap = BLUES.copy()
    cmap.set_bad("#e9e8e4")
    im = ax.imshow(shown, cmap=cmap, vmin=vmin, vmax=vmax)
    for i in range(n):
        for j in range(n):
            if i != j:
                dark = (mat[i][j] - vmin) / (vmax - vmin) > 0.6
                ax.text(j, i, f"{mat[i][j]:.2f}", ha="center", va="center",
                        fontsize=10, color="white" if dark else INK)
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("%s (scale starts at %.2f)" % (metric, vmin), fontsize=9, color=INK)
    cbar.outline.set_visible(False)
    ax.set_xticks(range(n), names, rotation=30, ha="right", fontsize=9, color=INK)
    ax.set_yticks(range(n), names, fontsize=9, color=INK)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    names, mat = compute_rf_matrix(trees)
    print_matrix(names, mat)
    save_csv(names, mat, args.out / "rf_matrix.csv")
    title = (f"Normalized Robinson-Foulds distance, {args.title} trees\n"
             "0 = identical topology, 1 = no shared splits (unrooted)")
    plot_heatmap(names, mat, args.out / "rf_heatmap.png", title)
    print("wrote", args.out / "rf_matrix.csv", "and", args.out / "rf_heatmap.png")


if __name__ == "__main__":
    main()
