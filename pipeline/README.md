# Generic PhyloBench pipeline (OB45, AG45, AS45, CH45, EB45)

One set of scripts for all five datasets; `X` is the two-letter dataset code (`OB`, `AG`, `AS`, `CH`, `EB`).
The original per-dataset scripts (`OB45/*.py`, `AG45/*.py`, `MsaPhylo/run_*45_*.py`, ...) are left as they were;
these are generalised copies and give identical files (checked on AS45_3, EB45_18 and OB45_7).
Commands are run from `DeepLearningClaudia/`.

## 1. Clean and concatenate

```bash
python pipeline/preprocess.py X                  # <X>45/*.afa      -> <X>45_cleaned/
python pipeline/concatenate_supermatrix.py X all # <X>45_cleaned/   -> <X>45_concat/<X>45_supermatrix.fasta
```

## 2. Find each method's boundary (largest number of loci it can still run on)

All subsets are random: `random.seed(42); random.sample(loci, N)`, as in `OB45/downsample_loci_7.py`.

| Method | Limit | Script |
|---|---|---|
| MSA Transformer (embedding + attention) | at most 1023 columns | `python pipeline/msat_boundary.py X` (computed) |
| CMAPLE | fails on too-divergent data ("Likelihood cost is very heavy") | `python pipeline/find_cmaple_boundary.py X [lo] [hi]` |
| Phyloformer | ~5.5 GB of RAM per 1000 columns (60 taxa), 1 TB machine ≈ 174k columns | `python pipeline/find_phyloformer_boundary.py X [lo] [hi]` |
| IQ-TREE | no hard limit, only run time | – |

Probe logs: `<X>45/logs/cmaple_boundary/`, `<X>45/logs/phyloformer_boundary/` (`probes.tsv`, `boundary.txt`).
With one seed, random.sample keeps the same order, so smaller samples are subsets of larger ones up to about N = 84 of 286-856 loci
and N = 340 of ~3000 loci, and only partly beyond that. CMAPLE can still succeed above a failure (a sample can fail while a larger
sample containing it succeeds); the boundary is the largest success directly below the smallest failure.  Phyloformer probes run one at a time and take
`Phyloformer/.big_run.lock`; other programs using memory at the same time lower its boundary.

## 3. Run the three stages

| Stage | Methods that run | N |
|---|---|---|
| 1 | all five | largest N that fits MSA-T (≤ 1023 columns) and on which CMAPLE also works |
| 2 | all but MSA-T | CMAPLE boundary |
| 3 | IQ-TREE + Phyloformer | Phyloformer boundary (or all loci if Phyloformer fits) |

```bash
pipeline/run_subset.sh X N [IQTREE_THREADS]        # downsample, concatenate, run all five methods
pipeline/collect_result_trees.sh X N                # -> <X>45_<N>_result_trees/
python pipeline/make_compare_R.py X N IQ,PF,CM,MA,ME   # methods that produced a tree
cd Comparison && Rscript compare_simple_<X>_<N>.R > figures_R/compare_simple_<X>_<N>.log 2>&1
                                                    # -> figures_R/compare_simple_<X>_<N>.pdf (+ RF/quartet numbers in the .log)
```

Method codes: IQ = IQ-TREE, PF = Phyloformer, CM = CMAPLE, MA = MSA-T attention tree, ME = MSA-T embedding tree.
Logs of each run: `<X>45/logs/subset_runs_<N>/`.

## Plots

`compare_template.R` is `Comparison/compare_simple_OB_603.R`'s analysis: every tree is rooted on the root tip of
`<X>-astral.tre` (`config.sh`), coloured by the NCBI clades in `clades/<X>.R`; RF and quartet distances use
unrooted trees.
