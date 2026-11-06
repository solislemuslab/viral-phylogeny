#!/usr/bin/env python3
"""Copy DeepLearningClaudia/ (minus the model folders) into a local clone of
github.com/solislemuslab/viral-phylogeny, keeping the same relative paths.

    python3 pipeline/stage_for_github.py [TARGET_REPO]     # default ../viral-phylogeny-github
    python3 pipeline/stage_for_github.py --dry-run         # list what would be copied

Left out:
  - the model / tool folders: CMAPLE, MsaPhylo, NeuralNJ, Phyloformer, Phyloformer2, PhyloGFN
  - the full GISAID download (gisaid_*.fasta, covid_all.fasta, covid_unique.fasta); the COVID subsets are kept
  - OB45.tar.gz (same content as OB45/), .DS_Store, __pycache__, .claude, lock files
  - any file > 50 MB (GitHub refuses files > 100 MB)

Files copied by a previous run that are no longer selected are removed from the target
(tracked in TARGET/.stage_manifest). Files of the GitHub repo that this script does not produce
(README.md, COVID_dataset/, PhyloBench/, Sample_data/, ...) are never touched.
"""
import fnmatch
import os
import shutil
import sys

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAX_BYTES = 50 * 1024 * 1024

MODEL_DIRS = {"CMAPLE", "MsaPhylo", "NeuralNJ", "Phyloformer", "Phyloformer2", "PhyloGFN"}
SKIP_DIR_NAMES = {".git", ".venv", "__pycache__", ".idea", ".claude", ".ipynb_checkpoints"}
SKIP_FILE_GLOBS = [".DS_Store", "*.pyc", "*.lock", "OB45.tar.gz",
                   "gisaid_*.fasta", "covid_all.fasta", "covid_unique.fasta"]


def selected():
    for root, dirs, files in os.walk(SRC):
        relroot = os.path.relpath(root, SRC)
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIR_NAMES
                         and not (relroot == "." and d in MODEL_DIRS))
        for f in sorted(files):
            if any(fnmatch.fnmatch(f, g) for g in SKIP_FILE_GLOBS):
                continue
            rel = f if relroot == "." else relroot + "/" + f
            p = os.path.join(SRC, rel)
            if os.path.islink(p):
                continue
            if os.path.getsize(p) > MAX_BYTES:
                print("  skip (>50 MB): " + rel, file=sys.stderr)
                continue
            yield rel


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry-run" in sys.argv
    dst = os.path.abspath(args[0] if args else os.path.join(SRC, "..", "viral-phylogeny-github"))
    files = list(selected())
    total = sum(os.path.getsize(os.path.join(SRC, r)) for r in files)
    if dry:
        print("\n".join(files))
        print("# %d files, %.1f MB" % (len(files), total / 1e6), file=sys.stderr)
        return
    if not os.path.isdir(os.path.join(dst, ".git")):
        sys.exit("%s is not a git clone; run: git clone "
                 "https://github.com/solislemuslab/viral-phylogeny.git %s" % (dst, dst))
    man = os.path.join(dst, ".stage_manifest")
    old = set(open(man).read().split("\n")) - {""} if os.path.exists(man) else set()
    new = set(files)
    for rel in sorted(old - new):
        p = os.path.join(dst, rel)
        if os.path.exists(p):
            os.remove(p)
    copied = 0
    for rel in files:
        s, d = os.path.join(SRC, rel), os.path.join(dst, rel)
        if os.path.exists(d) and os.path.getsize(d) == os.path.getsize(s) \
                and os.path.getmtime(d) >= os.path.getmtime(s):
            continue
        os.makedirs(os.path.dirname(d), exist_ok=True)
        shutil.copy2(s, d)
        copied += 1
    with open(man, "w") as fh:
        fh.write("\n".join(sorted(new)) + "\n")
    print("%d files (%.1f MB) selected, %d copied, %d stale removed -> %s"
          % (len(files), total / 1e6, copied, len(old - new), dst))


if __name__ == "__main__":
    main()
