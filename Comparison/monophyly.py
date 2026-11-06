"""
Are groups monophyletic? (item b.vii)

Groups can be taxonomic (as for the bundled PF00066 test trees:
PF00066_groups.tsv, broad group and subgroup per species) or, for a virus data
set, geographic region -- pass a file with `label<TAB>region[<TAB>subgroup]` rows:

    python monophyly.py --tree-dir /path/to/trees --groups regions.tsv

With no --groups (and none found next to --tree-dir), the group/subgroup levels
are skipped. Separately, when leaf labels look like Pfam-style
'ACCESSION_SPECIES/start-end' (tree_utils.has_domain_ranges), a "protein copies"
level checks whether multiple domain copies of the same accession stay together;
this is skipped when labels don't have that shape (e.g. viral accession ids).

A group is monophyletic in a tree when one branch separates all of its tips from
all other tips (unrooted, so the answer does not depend on the root).  When it is
not, "intruders" is the number of foreign tips inside the smallest clade that
still contains the whole group: 0 = monophyletic, larger = more mixed.

Writes to <out>/:
  monophyly.csv   group, size, tree, intruders, monophyletic (level = group / subgroup / protein)
  monophyly.png   groups x methods (+ tick = monophyletic, +k = k intruders) and how many
                  multi-copy proteins have their copies together (only the levels that apply)
"""
import csv

import matplotlib.pyplot as plt
import numpy as np

from compare import BLUES
from tree_utils import (INK, SURFACE, has_domain_ranges, load_all, load_groups, parse_args,
                        popcount, protein_of)


def intruders(tree, gmask):
    """Foreign tips in the smallest clade (either side of any split) holding the whole group."""
    full = (1 << tree.n) - 1
    best = tree.n
    for s in tree.splits:
        for side in (s, full ^ s):
            if side & gmask == gmask:
                best = min(best, popcount(side))
    return best - popcount(gmask)


def protein_masks(taxa):
    by = {}
    for i, t in enumerate(taxa):
        by[protein_of(t)] = by.get(protein_of(t), 0) | (1 << i)
    return {p: m for p, m in by.items() if popcount(m) >= 2}


def collect(trees, groups, taxa):
    """rows of (level, group, size, tree, intruders); only the levels that apply."""
    levels = []
    if groups is not None:
        levels += [("group", groups.masks("broad")), ("subgroup", groups.masks("sub"))]
    if has_domain_ranges(taxa):
        levels.append(("protein", protein_masks(taxa)))
    rows = []
    for level, masks in levels:
        for g, m in masks.items():
            for name, t in trees.items():
                rows.append((level, g, popcount(m), name, intruders(t, m)))
    return rows


def plot_group_panel(fig, rect, rows, names, groups, n):
    labelled = []
    for level in ("group", "subgroup"):
        gs = list(dict.fromkeys(r[1] for r in rows if r[0] == level))
        if level == "group":
            gs = [g for g in groups.order if g in gs]
        labelled += [(level, g) for g in gs]
    lookup = {(r[0], r[1], r[3]): r for r in rows}
    ax = fig.add_axes(rect, facecolor=SURFACE)
    grid = np.zeros((len(labelled), len(names)))
    for i, (level, g) in enumerate(labelled):
        for j, m in enumerate(names):
            _, _, size, _, k = lookup[(level, g, m)]
            grid[i, j] = k / (n - size) if n > size else 0
    ax.imshow(grid, cmap=BLUES, vmin=0, vmax=max(0.5, grid.max()), aspect="auto")
    for i, (level, g) in enumerate(labelled):
        for j, m in enumerate(names):
            k = lookup[(level, g, m)][4]
            dark = grid[i, j] / max(0.5, grid.max()) > 0.6
            ax.text(j, i, "✓" if k == 0 else f"+{k}", ha="center", va="center", fontsize=9,
                    color="white" if dark else INK, fontweight="bold" if k == 0 else "normal")
    ax.set_xticks(range(len(names)), names, rotation=30, ha="right", fontsize=9, color=INK)
    ax.set_yticks(range(len(labelled)),
                  [f"{g}  ({lookup[(lv, g, names[0])][2]})" for lv, g in labelled],
                  fontsize=9, color=INK)
    n_broad = sum(1 for lv, _ in labelled if lv == "group")
    if 0 < n_broad < len(labelled):
        ax.axhline(n_broad - 0.5, color=SURFACE, lw=5)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title("Monophyly of taxonomic groups\n✓ = monophyletic; +k = k foreign tips inside the "
                 "smallest clade holding the group;\n(group size in brackets; subgroups below the gap)",
                 loc="left", fontsize=9.5, color=INK)


def plot_protein_panel(fig, rect, rows, names, colors, title):
    prot = {m: sum(1 for r in rows if r[0] == "protein" and r[3] == m and r[4] == 0) for m in names}
    n_prot = len({r[1] for r in rows if r[0] == "protein"})
    ax = fig.add_axes(rect, facecolor=SURFACE)
    y = np.arange(len(names))[::-1]
    ax.barh(y, [prot[m] for m in names], color=[colors[m] for m in names], height=0.6)
    for yi, m in zip(y, names):
        ax.text(prot[m] + 0.4, yi, f"{prot[m]} of {n_prot}", va="center", fontsize=9, color=INK)
    ax.set_yticks(y, names, fontsize=9, color=INK)
    ax.set_xlim(0, n_prot * 1.25)
    ax.set_xticks([])
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title(f"Proteins with several {title} copies\nwhose copies form one clade",
                 loc="left", fontsize=9.5, color=INK)


def plot_grid(rows, names, groups, colors, title, n, path):
    has_group = any(r[0] in ("group", "subgroup") for r in rows)
    has_protein = any(r[0] == "protein" for r in rows)
    n_labelled = len({(r[0], r[1]) for r in rows if r[0] in ("group", "subgroup")})
    if has_group and has_protein:
        fig = plt.figure(figsize=(11.5, 0.42 * n_labelled + 2.6), facecolor=SURFACE)
        plot_group_panel(fig, [0.0, 0.0, 0.62, 1.0], rows, names, groups, n)
        plot_protein_panel(fig, [0.74, 0.0, 0.24, 1.0], rows, names, colors, title)
    elif has_group:
        fig = plt.figure(figsize=(8, 0.42 * n_labelled + 2.6), facecolor=SURFACE)
        plot_group_panel(fig, [0.0, 0.0, 0.98, 1.0], rows, names, groups, n)
    else:
        fig = plt.figure(figsize=(5.5, 3.6), facecolor=SURFACE)
        plot_protein_panel(fig, [0.0, 0.0, 0.98, 1.0], rows, names, colors, title)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main(args=None):
    args = args or parse_args(__doc__.strip().split("\n\n")[0])
    args.out.mkdir(parents=True, exist_ok=True)
    trees = load_all(args.tree_dir, args.reference)
    taxa = trees.taxa
    names = list(trees)
    groups = load_groups(taxa, args.groups, args.tree_dir, args.title)
    if groups is None:
        print("monophyly: no groups file found (see --groups); skipping the group/subgroup levels")
    if not has_domain_ranges(taxa):
        print("monophyly: labels don't look like Pfam-style 'ACCESSION_SPECIES/start-end'; "
              "skipping the multi-copy-protein level")
    rows = collect(trees, groups, taxa)
    if not rows:
        print("monophyly: nothing to check (no groups file and no domain-range labels); "
              "pass --groups to check group monophyly")
        return
    with open(args.out / "monophyly.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["level", "group", "size", "tree", "intruders", "monophyletic"])
        for r in rows:
            w.writerow(list(r) + ["yes" if r[4] == 0 else "no"])
    for level in ("group", "subgroup", "protein"):
        sel = [r for r in rows if r[0] == level]
        total = len({r[1] for r in sel})
        if not total:
            continue
        print(f"{level}: monophyletic groups out of {total}")
        for m in names:
            print(f"  {m:16s} {sum(1 for r in sel if r[3] == m and r[4] == 0)}")
    plot_grid(rows, names, groups, trees.colors, args.title, len(taxa), args.out / "monophyly.png")
    print("wrote monophyly.csv and monophyly.png in", args.out)


if __name__ == "__main__":
    main()
