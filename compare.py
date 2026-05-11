from ete3 import Tree
from Bio import Phylo
import matplotlib.pyplot as plt


def compute_rf_distance(tree1_path, tree2_path):
    """
    Compute RF distance using ETE3.
    """
    tree1 = Tree(tree1_path, format=1)
    tree2 = Tree(tree2_path, format=1)

    rf, max_rf, common_leaves, *_ = tree1.robinson_foulds(
        tree2,
        unrooted_trees=True
    )

    print("RF distance:", rf)
    print("Maximum RF distance:", max_rf)
    print("Normalized RF distance:", rf / max_rf if max_rf != 0 else 0)
    print("Number of common leaves:", len(common_leaves))


def plot_cophylogeny(tree1_path, tree2_path):
    """
    Draw a simple co-phylogeny
    """
    # Load trees with Biopython, not ETE3
    tree1 = Phylo.read(tree1_path, "newick")
    tree2 = Phylo.read(tree2_path, "newick")

    # Get shared leaf names only
    tips1 = {t.name for t in tree1.get_terminals()}
    tips2 = {t.name for t in tree2.get_terminals()}
    common = sorted(tips1 & tips2)

    print("Common taxa:", len(common))

    # Prune taxa not shared
    for tip in list(tree1.get_terminals()):
        if tip.name not in common:
            tree1.prune(tip)

    for tip in list(tree2.get_terminals()):
        if tip.name not in common:
            tree2.prune(tip)

    # Visualize
    y1 = {tip.name: i for i, tip in enumerate(tree1.get_terminals())}
    y2 = {tip.name: i for i, tip in enumerate(tree2.get_terminals())}

    fig = plt.figure(figsize=(14, 20))
    ax1 = fig.add_axes([0.02, 0.05, 0.35, 0.9]) # Left region
    ax2 = fig.add_axes([0.63, 0.05, 0.35, 0.9]) # right region
    axm = fig.add_axes([0.37, 0.05, 0.26, 0.9]) # middle region

    Phylo.draw(tree1, do_show=False, axes=ax1, show_confidence=False)
    Phylo.draw(tree2, do_show=False, axes=ax2, show_confidence=False)

    ax1.set_title("Tree 1")
    ax2.set_title("Tree 2")
    axm.set_xlim(0, 1) # Set left and right edge for middle region (used for later connected lines)
    axm.set_ylim(-1, len(common))
    axm.axis("off")

    # Draw connecting lines between identical taxa
    for name in common:
        axm.plot([0, 1], [y1[name], y2[name]], lw=0.5)

    plt.show()


if __name__ == "__main__":
    tree1_path = "../IQTREE_results/iqtree.nwk"
    tree2_path = "../CMAPLE/CMAPLE_results/cmaple.nwk"

    compute_rf_distance(tree1_path, tree2_path)
    plot_cophylogeny(tree1_path, tree2_path)
