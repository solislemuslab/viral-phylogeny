"""
Shared helpers for the tree comparison scripts. Nothing here is specific to one
Pfam family or tree set -- every script takes `--tree-dir` (and optionally
`--reference`, `--groups`, `--out`, `--title`; see `parse_args`/`build_arg_parser`
below) and works the same way for any directory of trees over the same taxa.

Locating the trees (`discover_trees`, used by `load_all`):
  - if `<tree-dir>/trees.tsv` (or `manifest.tsv`) exists, it is read as a
    tab-separated table with columns `file`, `method` (display name) and an
    optional `reference` column ("yes" on one row); this is how the bundled
    PF00066_result_trees/trees.tsv names its five trees.
  - otherwise every *.nwk/*.newick/*.treefile/*.tree/*.nex file in the directory
    is used, named after its filename with any prefix shared by all of them
    stripped (e.g. "PF00066_IQTree.treefile" -> "IQTree").
  - the reference tree (used for midpoint rooting and as the fixed point every
    other tree is compared to) is `--reference`, else the manifest's, else
    whichever tree sorts first.

Every tree is reduced to its set of non-trivial splits (bipartitions stored as
integer bitmasks over a common, sorted taxon list -- the reference tree's
leaves; every other tree must have exactly the same leaf set).  That is all RF
distance, quartets, consensus, reference-clade and monophyly checks need, and it
makes polytomies explicit: an edge whose length is <= COLLAPSE_TOL is dropped,
so an unresolved node simply has fewer splits than a binary one.

For drawing, a rooted view is rebuilt from the splits plus a root split
(rooted_tree).  The reference tree is midpoint-rooted; every other tree is
rooted on the edge whose split is most similar to the reference root, so the
trees are drawn in a comparable orientation.  Drawings are cladograms (branch
lengths are not comparable across methods: ML substitutions, embedding
distances, placeholder lengths).

Taxon groups (`load_groups`, item b.vii and the group colouring in the other
plots) are optional: a tab-separated `key, group[, subgroup]` file, found via
`--groups`, or `groups.tsv`/`taxon_groups.tsv` in the tree directory, or
`<title>_groups.tsv` next to this file (e.g. PF00066_groups.tsv is the bundled
example).  `key` may be a full leaf label or, for Pfam-style
`ACCESSION_SPECIES/start-end` labels, the species code -- see `species_of`.
When no groups file is found, `load_groups` returns None and every script
degrades to plain (ungrouped) styling instead of failing.
"""
import argparse
import csv
import os
import re
from pathlib import Path

import dendropy

HERE = Path(__file__).resolve().parent
# only a CLI default so the scripts run out of the box on the bundled test data;
# nothing below depends on this particular directory.
DEFAULT_TREE_DIR = HERE.parent / "PF00066_result _trees"

COLLAPSE_TOL = 1e-6  # internal edges this short (or shorter) count as polytomies
TREE_EXTENSIONS = {".nwk", ".newick", ".treefile", ".tree", ".nex", ".nexus"}
MANIFEST_NAMES = ("trees.tsv", "manifest.tsv")

# ---- palette (dataviz reference palette, light surface) -------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8983"
NEUTRAL = "#c3c2b7"   # "nothing to report" edges and cells
GRID = "#e9e8e4"
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
               "#008300", "#4a3aa7", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]  # secondary encoding for groups
_BLUE_STOPS = [(0.0, "#cde2fb"), (0.5, "#3987e5"), (1.0, "#0d366b")]  # sequential ramp


def popcount(x):
    return bin(x).count("1")


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def _rgb_to_hex(rgb):
    return "#%02x%02x%02x" % tuple(max(0, min(255, round(c))) for c in rgb)


def blues(k):
    """k ordinal steps, light -> dark, along the sequential blue ramp (k >= 1)."""
    if k <= 1:
        return ["#2a78d6"]
    out = []
    for i in range(k):
        t = i / (k - 1)
        (t0, c0), (t1, c1) = next((a, b) for a, b in zip(_BLUE_STOPS, _BLUE_STOPS[1:])
                                  if a[0] <= t <= b[0])
        f = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
        r0, r1 = _hex_to_rgb(c0), _hex_to_rgb(c1)
        out.append(_rgb_to_hex(tuple(a + (b - a) * f for a, b in zip(r0, r1))))
    return out


def method_colors(names):
    """{tree name: colour}, fixed categorical order; MUTED past the 8th."""
    return {n: (CATEGORICAL[i] if i < len(CATEGORICAL) else MUTED) for i, n in enumerate(names)}


# --------------------------------------------------------------------------
# labels
# --------------------------------------------------------------------------
def species_of(label):
    """Pfam-style 'ACCESSION_SPECIES/start-end' -> 'SPECIES'.

    A label without an '_' (e.g. a viral strain/accession id) is returned
    unchanged, so it can be matched directly against a groups file keyed by
    the full label.
    """
    head = label.split("/", 1)[0]
    return head.split("_", 1)[1] if "_" in head else head


def protein_of(label):
    """The part of a label before an optional '/start-end' domain-range suffix.

    Several rows can share this (Pfam domain copies of one protein); a label
    with no '/' is its own protein_of.
    """
    return label.split("/", 1)[0]


def has_domain_ranges(taxa):
    """True if labels look like Pfam-style 'ACCESSION_SPECIES/start-end'."""
    return any("/" in t for t in taxa)


# --------------------------------------------------------------------------
# locating the trees
# --------------------------------------------------------------------------
def _read_manifest(path):
    files, reference = {}, None
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            fname = (row.get("file") or "").strip()
            if not fname:
                continue
            name = (row.get("method") or Path(fname).stem).strip()
            files[name] = fname
            if (row.get("reference") or "").strip().lower() in ("yes", "true", "1"):
                reference = name
    if not files:
        raise ValueError(f"manifest {path} has no rows (need a 'file' column)")
    return files, reference


def _auto_names(stems):
    """Display name per stem: the shared leading '..._' token stripped off."""
    prefix = os.path.commonprefix(stems)
    cut = max(prefix.rfind("_"), prefix.rfind("-"))
    prefix = prefix[:cut + 1] if cut >= 0 else ""
    names, seen = {}, set()
    for s in stems:
        name = s[len(prefix):].strip("_- ") or s
        if name in seen:  # keep names unique even if stripping collided two stems
            name = s
        seen.add(name)
        names[s] = name
    return names


def discover_trees(tree_dir):
    """({display name: filename}, reference name or None), in display order.

    See the module docstring for the manifest / auto-discovery rules.
    """
    tree_dir = Path(tree_dir)
    for mname in MANIFEST_NAMES:
        mpath = tree_dir / mname
        if mpath.exists():
            return _read_manifest(mpath)
    paths = sorted(p for p in tree_dir.iterdir()
                   if p.is_file() and p.suffix.lower() in TREE_EXTENSIONS)
    if not paths:
        raise FileNotFoundError(
            f"no tree files ({', '.join(sorted(TREE_EXTENSIONS))}) in {tree_dir}; "
            "add one, or a trees.tsv manifest with columns file, method[, reference]")
    auto = _auto_names([p.stem for p in paths])
    return {auto[p.stem]: p.name for p in paths}, None


def dataset_title(tree_dir):
    """A short dataset name for figure titles, derived from the directory name."""
    name = Path(tree_dir).name.strip()
    name = re.sub(r"[ _]*result[ _]*trees$", "", name, flags=re.I)
    name = re.sub(r"[ _]*trees$", "", name, flags=re.I)
    return name.strip(" _") or Path(tree_dir).name


def default_fig_dir(tree_dir):
    """Comparison/figures/<tree-dir name>, so different tree sets don't collide."""
    slug = re.sub(r"\s+", "_", Path(tree_dir).name.strip())
    return HERE / "figures" / slug


# --------------------------------------------------------------------------
# CLI plumbing shared by every script
# --------------------------------------------------------------------------
def build_arg_parser(description=None):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--tree-dir", type=Path, default=DEFAULT_TREE_DIR,
                   help="directory of tree files (.nwk/.newick/.treefile/.tree/.nex), "
                        "optionally with a trees.tsv manifest (file, method[, reference]); "
                        f"default: {DEFAULT_TREE_DIR.name} (the bundled test trees)")
    p.add_argument("--reference", default=None,
                   help="tree to use as the reference (default: the manifest's, "
                        "else whichever tree sorts first)")
    p.add_argument("--groups", type=Path, default=None,
                   help="tab-separated taxon-group file (key, group[, subgroup]); "
                        "default: groups.tsv/taxon_groups.tsv in --tree-dir, or "
                        "<title>_groups.tsv next to tree_utils.py, if present")
    p.add_argument("--out", type=Path, default=None,
                   help="output directory for figures/tables (default: "
                        "figures/<tree-dir name> next to tree_utils.py)")
    p.add_argument("--title", default=None,
                   help="dataset name used in figure titles (default: the tree-dir "
                        "name, with a trailing '_trees'/'_result_trees' stripped)")
    return p


def parse_args(description=None, argv=None):
    args = build_arg_parser(description).parse_args(argv)
    args.tree_dir = Path(args.tree_dir)
    args.title = args.title or dataset_title(args.tree_dir)
    args.out = Path(args.out) if args.out else default_fig_dir(args.tree_dir)
    return args


# --------------------------------------------------------------------------
# unrooted trees as split sets
# --------------------------------------------------------------------------
class Tree:
    """An unrooted topology: `splits` maps a normalized bitmask -> (length, label).

    A mask is normalized to the side that does NOT contain taxon 0.
    """

    def __init__(self, name, taxa, splits, n_collapsed, midpoint_mask=None):
        self.name = name
        self.taxa = taxa
        self.n = len(taxa)
        self.splits = splits
        self.n_collapsed = n_collapsed
        self.midpoint_mask = midpoint_mask

    @property
    def resolution(self):
        """Fraction of the n-3 splits a fully binary tree would have."""
        return len(self.splits) / (self.n - 3)


class TreeSet(dict):
    """{method name: Tree}, plus the dataset-wide info the scripts share."""
    reference = None   # resolved reference tree name
    taxa = None         # the common, sorted taxon list
    tree_dir = None
    colors = None        # {name: colour}, fixed categorical order


def _leaf_masks(tree, index):
    mask = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            mask[nd] = 1 << index[nd.taxon.label]
        else:
            m = 0
            for c in nd.child_nodes():
                m |= mask[c]
            mask[nd] = m
    return mask


def _normalize(m, full):
    return full ^ m if m & 1 else m


def _midpoint_mask(path, index):
    tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    tree.reroot_at_midpoint(update_bipartitions=False)
    full = (1 << len(index)) - 1
    mask = _leaf_masks(tree, index)
    return _normalize(mask[tree.seed_node.child_nodes()[0]], full)


def _parse(name, path, taxa, index, tol):
    tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    these = {nd.taxon.label for nd in tree.leaf_node_iter()}
    if these != set(index):
        missing, extra = set(index) - these, these - set(index)
        detail = []
        if missing:
            detail.append(f"missing {len(missing)} (e.g. {sorted(missing)[:3]})")
        if extra:
            detail.append(f"has {len(extra)} extra (e.g. {sorted(extra)[:3]})")
        raise ValueError(f"{name} ({path.name}) has different taxa than the reference: "
                         + "; ".join(detail))
    n = len(taxa)
    full = (1 << n) - 1
    mask = _leaf_masks(tree, index)
    root = tree.seed_node
    nonroot = [nd for nd in tree.preorder_node_iter() if nd is not root]
    has_lengths = any(nd.edge_length is not None for nd in nonroot)
    edges = {}
    for nd in nonroot:
        if nd.is_leaf():
            continue
        s = _normalize(mask[nd], full)
        if popcount(s) < 2 or n - popcount(s) < 2:
            continue
        length = (nd.edge_length or 0.0) if has_lengths else 1.0
        label = nd.label.strip() if nd.label else None
        if s in edges:  # rooted bifurcating root: two edges, one split
            length += edges[s][0]
            label = edges[s][1] or label
        edges[s] = (length, label)
    splits = {s: e for s, e in edges.items() if e[0] > tol}
    return Tree(name, taxa, splits, len(edges) - len(splits))


def load_all(tree_dir, reference=None, tol=COLLAPSE_TOL):
    """TreeSet: {method name: Tree}, over the tree files found in `tree_dir`.

    Every tree must share exactly the reference tree's taxa (see `_parse`).
    `reference` overrides the manifest's / the first tree's choice; it is
    midpoint-rooted, and every other tree is rooted to match it (root_masks).
    """
    tree_dir = Path(tree_dir)
    files, manifest_ref = discover_trees(tree_dir)
    reference = reference or manifest_ref or next(iter(files))
    if reference not in files:
        raise KeyError(f"--reference {reference!r} is not one of the trees: "
                       f"{', '.join(files)}")
    ref_path = tree_dir / files[reference]
    ref_dp = dendropy.Tree.get(path=str(ref_path), schema="newick", preserve_underscores=True)
    taxa = sorted(nd.taxon.label for nd in ref_dp.leaf_node_iter())
    index = {t: i for i, t in enumerate(taxa)}
    trees = TreeSet()
    for name, fname in files.items():
        trees[name] = _parse(name, tree_dir / fname, taxa, index, tol)
    trees[reference].midpoint_mask = _midpoint_mask(ref_path, index)
    trees.reference, trees.taxa, trees.tree_dir = reference, taxa, tree_dir
    trees.colors = method_colors(list(files))
    return trees


def rf_distance(a, b):
    """(RF, maximum possible RF) between two unrooted trees.

    max = |splits(a)| + |splits(b)|, so a tree with polytomies is not
    penalised against the 2(n-3) of a fully binary pair.
    """
    return len(a.splits.keys() ^ b.splits.keys()), len(a.splits) + len(b.splits)


def compatible(a, b):
    """Two normalized splits can co-exist in one tree."""
    return not (a & b) or not (a & ~b) or not (b & ~a)


# --------------------------------------------------------------------------
# rooting and the rooted view
# --------------------------------------------------------------------------
def best_match_root(masks, ref_root, n):
    """The split in `masks` most similar to `ref_root` (None if there is none)."""
    full = (1 << n) - 1
    rc = full ^ ref_root
    cands = set(masks)
    if popcount(ref_root) == 1 or popcount(rc) == 1:
        cands.add(ref_root)  # a leaf edge is a valid root only if the reference has it
    best = None
    for s in cands:
        score = min(popcount(s ^ ref_root), popcount(s ^ rc))
        balance = min(popcount(s), n - popcount(s))
        key = (score, -balance, s)
        if best is None or key < best[0]:
            best = (key, s)
    return best[1] if best else None


def root_masks(trees, reference=None):
    """{name: root split} - reference midpoint-rooted, others matched to it."""
    reference = reference or trees.reference
    ref = trees[reference]
    out = {reference: ref.midpoint_mask}
    for name, t in trees.items():
        if name != reference:
            out[name] = best_match_root(t.splits, ref.midpoint_mask, t.n)
    return out


class RNode:
    __slots__ = ("mask", "children", "leaf", "height", "size")

    def __init__(self, mask):
        self.mask = mask
        self.children = []
        self.leaf = None      # taxon index for a leaf
        self.height = 0       # edges to the deepest leaf below (leaf = 0)
        self.size = 1         # leaves below

    def split(self, n):
        return _normalize(self.mask, (1 << n) - 1)

    def leaves(self):
        if self.leaf is not None:
            return [self]
        return [l for c in self.children for l in c.leaves()]

    def walk(self):
        yield self
        for c in self.children:
            yield from c.walk()


def rooted_tree(splits, n, root_mask=None):
    """Rebuild a rooted tree (RNode) from unrooted splits and a root split.

    Each edge's clade is the side of its split lying away from the root edge.
    root_mask=None roots on taxon 0's leaf edge (only a fallback).
    """
    full = (1 << n) - 1
    clades = {1 << i for i in range(n)}
    clades.add(full)
    if root_mask is None:
        clades.update(splits)
    else:
        clades.update((root_mask, full ^ root_mask))
        for s in splits:
            if s == root_mask:
                continue
            for x in (s, full ^ s):
                if not (x & ~root_mask) or not (x & root_mask):
                    clades.add(x)
                    break
            else:
                raise ValueError("split incompatible with the root split")
    order = sorted(clades, key=lambda m: (popcount(m), m))
    nodes = {m: RNode(m) for m in order}
    for i, m in enumerate(order[:-1]):
        for p in order[i + 1:]:
            if m & p == m:
                nodes[p].children.append(nodes[m])
                break
    for m, nd in nodes.items():
        if popcount(m) == 1:
            nd.leaf = m.bit_length() - 1
    root = nodes[full]
    _annotate(root)
    return root


def _annotate(nd):
    if nd.children:
        for c in nd.children:
            _annotate(c)
        nd.height = 1 + max(c.height for c in nd.children)
        nd.size = sum(c.size for c in nd.children)


def ladderize(nd):
    """Sort children small clade first (ties: lowest taxon index); in place."""
    for c in nd.children:
        ladderize(c)
    nd.children.sort(key=lambda c: (c.size, c.mask & -c.mask))
    return nd


def to_newick(nd, taxa, label_fn=None):
    if nd.leaf is not None:
        return taxa[nd.leaf]
    inner = ",".join(to_newick(c, taxa, label_fn) for c in nd.children)
    lab = label_fn(nd) if label_fn and nd.mask != (1 << len(taxa)) - 1 else ""
    return "(" + inner + ")" + (str(lab) if lab is not None else "")


def rooted_all(trees, reference=None):
    """{name: ladderized RNode} for every tree, with the matched roots."""
    roots = root_masks(trees, reference)
    return {name: ladderize(rooted_tree(t.splits, t.n, roots[name])) for name, t in trees.items()}


# --------------------------------------------------------------------------
# taxon groups (taxonomy for PF00066; geographic region for a virus data set)
# --------------------------------------------------------------------------
class Groups:
    """label -> group (and optional subgroup) plus a fixed colour/marker per group."""

    def __init__(self, taxa, broad, sub, order):
        self.taxa = taxa
        self.broad = broad            # label -> group
        self.sub = sub                # label -> subgroup ('' if none)
        self.order = order            # group order = colour slot order

    def color(self, group):
        i = self.order.index(group)
        return CATEGORICAL[i] if i < len(CATEGORICAL) else MUTED

    def marker(self, group):
        i = self.order.index(group)
        return MARKERS[i] if i < len(MARKERS) else "o"

    def masks(self, level="broad"):
        """{group: bitmask} for groups with >= 2 members."""
        table = self.broad if level == "broad" else self.sub
        out = {}
        for i, t in enumerate(self.taxa):
            g = table[t]
            if g:
                out[g] = out.get(g, 0) | (1 << i)
        return {g: m for g, m in out.items() if popcount(m) >= 2}

    def group_of_mask(self, mask):
        """The single broad group that contains every leaf of `mask`, else None."""
        for g, m in self.masks("broad").items():
            if mask & ~m == 0:
                return g
        return None


def find_groups_file(tree_dir, title=None):
    """The first of the conventional groups-file locations that exists, else None."""
    tree_dir = Path(tree_dir)
    title = title or dataset_title(tree_dir)
    for c in (tree_dir / "groups.tsv", tree_dir / "taxon_groups.tsv",
             HERE / f"{title}_groups.tsv"):
        if c.exists():
            return c
    return None


def load_groups(taxa, path=None, tree_dir=None, title=None, required=False):
    """Read a tab-separated file: key, group[, subgroup] -> Groups, or None.

    The key is a full leaf label or a species code (see `species_of`); for a
    virus data set give leaf label -> region.  `path` wins; otherwise it is
    looked up from `tree_dir`/`title` via `find_groups_file`.  When no file is
    found: None if `required` is false (the caller draws ungrouped), else
    FileNotFoundError.  A taxon missing from a *found* file is treated as
    ungrouped, with a warning -- it does not block the rest.
    """
    if path is None and tree_dir is not None:
        path = find_groups_file(tree_dir, title)
    if path is None:
        if required:
            raise FileNotFoundError(
                "no taxon-group file found; pass --groups, or add groups.tsv / "
                "taxon_groups.tsv to the tree directory")
        return None
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"--groups file not found: {path}")
    rows, order = {}, []
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            rows[parts[0]] = (parts[1], parts[2] if len(parts) > 2 else "")
            if parts[1] not in order:
                order.append(parts[1])
    broad, sub, missing = {}, {}, []
    for t in taxa:
        hit = rows.get(t) or rows.get(species_of(t))
        if hit is None:
            missing.append(t)
            broad[t], sub[t] = "", ""
        else:
            broad[t], sub[t] = hit
    if missing:
        print(f"load_groups: {path.name} has no group for {len(missing)} of {len(taxa)} taxa "
              f"(e.g. {', '.join(missing[:5])}{' ...' if len(missing) > 5 else ''}); "
              "treating them as ungrouped")
    return Groups(taxa, broad, sub, order)
