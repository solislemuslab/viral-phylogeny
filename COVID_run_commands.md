# COVID 300: commands run to build and compare the trees

Run on 2026-10-22, CPU-only machine (no GPU), from the project root
`/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia` unless a `cd` says otherwise.

Input for every method: `COVID_data_fasta/covid_300_aligned.fasta`
(300 stratified SARS-CoV-2 genomes + 4 animal-CoV outgroups, MAFFT-aligned; 304 taxa x 30,032 sites).
That file was made earlier with the README preprocessing steps (cat, `remove_duplicates.py`,
`select_covid_300_stratified.py`, add outgroups, `mafft --auto`); those steps were not rerun.

Final trees (trees only, no logs): `covid_result_trees/`

| Method   | Tree in `covid_result_trees/` | Raw output and log |
|----------|-------------------------------|--------------------|
| IQ-TREE  | `covid_300_IQTree.nwk`        | `IQTREE_results/covid_300_aligned.fasta.*` |
| CMAPLE   | `covid_300_CMAPLE.treefile`   | `CMAPLE/CMAPLE_results/covid_300_aligned.fasta.*` |
| NeuralNJ | `covid_300_NeuralNJ.nwk`      | `NeuralNJ/NeuralNJ_results/` |

## Differences from the README

| Step | README | Run here | Why |
|------|--------|----------|-----|
| IQ-TREE | `-nt AUTO`, no seed, run inside `COVID_data_fasta/` then `mv` | `-T 8 --seed 603 --prefix IQTREE_results/...` | 8 threads to share the machine with other jobs; `--prefix` writes straight to the results folder. The seed is arbitrary (603 has no meaning for COVID). Model is the same (GTR+G). |
| CMAPLE | `-aln` only | adds `-st DNA --seed 603 --prefix ...` | same reasons; default model is GTR either way |
| NeuralNJ env | `conda env create -f environment.yaml` (CUDA) | new CPU env `NeuralNJ` | no GPU here |
| raxmlpy | `python setup.py build/install` | compiled the binding by hand against raxml-ng 1.2.2 | the shipped `libraxml.so` needs glibc 2.38 (built elsewhere) and the bundled raxml-ng source is 2.0 (coraxlib), which `raxmlpy.cpp` does not compile against |
| NeuralNJ input | `convert_to_phylip.py` output in `examples/` | `prepare_covid_phy.py` output in `examples/covid_300_aligned/` | NeuralNJ needs taxa named `T1..T304` and only encodes `ACGT-N` (104 IUPAC codes set to `N`); the config's `instance_path` is a folder |
| NeuralNJ run | `python finetune_rl_search.py ...` | same arguments, launched through `run_neuralnj_cpu.py` | the checkpoint was saved on CUDA; the wrapper loads it on CPU. `finetune_rl_search.py` is unchanged |

The old COVID outputs in `IQTREE_results/` (no `.treefile`, only a renamed `iqtree.nwk`) and
`CMAPLE/CMAPLE_results/` (log was an error: "treefile already exists") were deleted before rerunning.

## 1. Environments

```bash
CONDA=/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/Miniforge3-26.1.1-3/bin/conda

# NeuralNJ (CPU)
$CONDA create -y -n NeuralNJ -c conda-forge python=3.9 pybind11
~/.conda/envs/NeuralNJ/bin/pip install --index-url https://download.pytorch.org/whl/cpu torch==2.2.2
~/.conda/envs/NeuralNJ/bin/pip install "numpy<2" einops dendropy ete3 biopython tqdm pyyaml fvcore tensorboard six

# R for the comparison script
$CONDA create -y -n Rphylo -c conda-forge r-base=4.4 r-ape r-phangorn r-phytools r-mass
```

IQ-TREE uses the existing `iq_tree` env (`~/.conda/envs/iq_tree/bin/iqtree`);
CMAPLE uses the existing binary `CMAPLE/cmaple-1.1.0-Linux-intel/bin/cmaple`.

### raxmlpy for NeuralNJ

```bash
# raxml-ng 1.2.2 (pll-modules based) built as a shared library, outside the project
git clone NeuralNJ/RAxMLpy/raxml-ng raxml-ng-1.2.2
cd raxml-ng-1.2.2
git checkout 1.2.2
git submodule update --init --recursive
mkdir build && cd build
cmake .. -DBUILD_AS_LIBRARY=ON -DUSE_TERRAPHAST=OFF
make -j16
cd ../..

# install the library and the pybind11 binding into the NeuralNJ env
E=~/.conda/envs/NeuralNJ
PY=$E/bin/python
SP=$($PY -c "import site; print(site.getsitepackages()[0])")
SUF=$($PY -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
RP=NeuralNJ/RAxMLpy
cp raxml-ng-1.2.2/build/lib/raxml-ng.so $E/lib/libraxml.so
mkdir -p $SP/raxmlpy
cp $RP/raxmlpy/__init__.py $RP/raxmlpy/core.py $SP/raxmlpy/
g++ -O3 -shared -fPIC -std=c++14 -w $($PY -m pybind11 --includes) \
    -Iraxml-ng-1.2.2/src -Iraxml-ng-1.2.2/build/localdeps/include \
    $RP/cpp/raxmlpy.cpp -L$E/lib -Wl,-rpath,$E/lib -lraxml -lm \
    -o $SP/raxmlpy/cpp_binding$SUF
```

Check: NeuralNJ on its own example `examples/len1024taxa50/G_l_1024_n_50_0_0.02_71.phy`
gave normalised RF = 0 against the true tree.

## 2. IQ-TREE

```bash
~/.conda/envs/iq_tree/bin/iqtree -s COVID_data_fasta/covid_300_aligned.fasta \
    -m GTR+G -T 8 --seed 603 --prefix IQTREE_results/covid_300_aligned.fasta
cp IQTREE_results/covid_300_aligned.fasta.treefile covid_result_trees/covid_300_IQTree.nwk
```

## 3. CMAPLE

```bash
CMAPLE/cmaple-1.1.0-Linux-intel/bin/cmaple -aln COVID_data_fasta/covid_300_aligned.fasta \
    -st DNA --prefix CMAPLE/CMAPLE_results/covid_300_aligned.fasta --seed 603
cp CMAPLE/CMAPLE_results/covid_300_aligned.fasta.treefile covid_result_trees/covid_300_CMAPLE.treefile
```

CMAPLE warns that the sequences look too divergent for it (because of the outgroups) and recommends IQ-TREE; it still finishes (~8 s).

## 4. NeuralNJ

`config/finetune_reinforce_search_example.yaml` already has `instance_path: "./examples/covid_300_aligned/"`.

```bash
cd NeuralNJ
~/.conda/envs/NeuralNJ/bin/python prepare_covid_phy.py   # writes examples/covid_300_aligned/{covid_300_aligned.phy,taxon_map.tsv}
~/.conda/envs/NeuralNJ/bin/python run_neuralnj_cpu.py finetune_rl_search.py \
    --config ./config/finetune_reinforce_search_example.yaml --infer_opt Argmax \
    > covid_300_NeuralNJ.log 2>&1
# output: output/Argmax_dim64_patch1/covid_300_aligned.tre (tips are T1..T304)
mkdir -p NeuralNJ_results
mv covid_300_NeuralNJ.log NeuralNJ_results/
cp output/Argmax_dim64_patch1/covid_300_aligned.tre NeuralNJ_results/
cd ..
cp NeuralNJ/NeuralNJ_results/covid_300_aligned.tre covid_result_trees/covid_300_NeuralNJ.nwk
```

The NeuralNJ tree keeps the `T1..T304` names; the comparison script maps them back with
`NeuralNJ/examples/covid_300_aligned/taxon_map.tsv`. On CPU it needs >100 GB RAM for 304 taxa.

**Status (2026-10-22): no NeuralNJ tree yet.** The CPU run was stopped by hand after 19.6 min without
finishing (the progress bar stays at 0/1 until the whole tree is built). It used 60-118 GB RAM
while the EB45 Phyloformer boundary run (`infer_alns.py`, N642) was using ~808 GB, which left
<60 GB free; it was stopped so an out-of-memory kill would not hit the boundary run and record a
false boundary. NeuralNJ cannot resume, so it was restarted from scratch with a watcher that kills
only NeuralNJ if available memory drops below 40 GB. `solislemus-001` has no GPU.
That restart was also stopped (after 1 h) once the memory measurements in section 6 showed that
304 taxa x 30032 sites needs ~4 TB: **NeuralNJ cannot run covid_300 on this machine.** The log of
that run is kept in `NeuralNJ/NeuralNJ_results/covid_300_NeuralNJ.log` (stuck at 0/1, `exit=143`
after 3651 s); the first run's log was removed. There is no NeuralNJ tree for covid_300.
See section 6 for covid_136.

## 5. Comparison

`Comparison/compare_simple_COVID_300.R` is a new script adapted from `compare_simple_OB_603.R` (which is unchanged):

- no reference trees: the methods (IQ-TREE, CMAPLE, NeuralNJ) are compared with each other, and IQ-TREE is the reference
- rooted on the bat virus RaTG13 (`EPI_ISL_402131`); RF is computed unrooted
- branches shorter than 1e-5 are collapsed (one substitution = 1/30032 = 3.3e-5; IQ-TREE reports unsupported branches as 1e-6)
- tips coloured by continent (`covid_300_stratified_metadata.csv`), outgroups in red, labels shortened to the EPI_ISL accession
- quartets: a random sample of 200,000 (seed 42) instead of all choose(304, 4) ~ 3.5e8
- circular trees drawn without branch lengths (outgroup branches are ~1000x longer); a phylogram without the outgroups shows the lengths inside SARS-CoV-2
- PCoA needs ≥ 3 trees and NMDS needs ≥ 4 trees; with fewer they are skipped

```bash
cd Comparison
~/.conda/envs/Rphylo/bin/Rscript compare_simple_COVID_300.R
# -> figures_R/compare_simple_COVID_300.pdf, figures_R/quartets_common_COVID_300.csv
```

The script compares whichever of the three trees exist. Results with IQ-TREE and CMAPLE (2026-10-22):

| | IQ-TREE vs CMAPLE |
|---|---|
| normalised RF (unrooted) | 0.172 |
| quartet distance (200,000 sampled quartets) | 0.001 |
| sampled quartets resolved identically by both | 69.3% |
| quartets resolved (rest are in polytomies) | IQ-TREE 71.2%, CMAPLE 70.4% |
| IQ-TREE clades recovered by CMAPLE | 83.1% |

IQ-TREE log-likelihood: -77268.725 (4 min 37 s).

## 6. NeuralNJ limit and the downsampled covid_136

### Why covid_300 does not fit

At inference NeuralNJ does not compress the alignment: `load_pi_instance` calls `only_padding_sample`,
which keeps every column (the 1024-column compression `compress_pad_seq` is only used in training).
In the first decoding step it builds a (columns x 64) tensor for every pair of taxa, several times over
(`aggregate` in `model.py`), so peak memory grows with taxa^2 x columns.

Measured on covid_300 with taxa and columns sampled at random (seed 42), CPU, Argmax:

| taxa | columns | time | peak RAM |
|------|---------|------|----------|
| 304 | 1024 | 720 s | 141.9 GB |
| 304 | 2048 | 1253 s | 281.9 GB |
| 152 | 1024 | 333 s | 35.9 GB |

Peak RAM ~ 141.9 GB x (columns / 1024) x (taxa / 304)^2, plus ~58 GB for the encoder's attention
over 30,000 sites (estimated from the code, 2 x 8 heads x columns^2 x 4 bytes; not measured).
304 taxa x 30032 columns would need ~4.2 TB. With a 900 GB budget (machine has 1007 GB) and full
genomes, the limit is **136 taxa** (predicted peak ~889 GB).

### covid_136

4 outgroups + 132 of the 300 SARS-CoV-2 sequences (`random.seed(42); random.sample`), all-gap
columns removed: 136 taxa x 29,994 columns. covid_300 files and trees are unchanged.

```bash
cd COVID_data_fasta
python downsample_covid.py 136          # -> covid_136_aligned.fasta
cd ..

# IQ-TREE and CMAPLE: same settings as covid_300
CMAPLE/cmaple-1.1.0-Linux-intel/bin/cmaple -aln COVID_data_fasta/covid_136_aligned.fasta \
    -st DNA --prefix CMAPLE/CMAPLE_results/covid_136_aligned.fasta --seed 603
~/.conda/envs/iq_tree/bin/iqtree -s COVID_data_fasta/covid_136_aligned.fasta \
    -m GTR+G -T 8 --seed 603 --prefix IQTREE_results/covid_136_aligned.fasta
mkdir -p covid_136_result_trees
cp CMAPLE/CMAPLE_results/covid_136_aligned.fasta.treefile covid_136_result_trees/covid_136_CMAPLE.treefile
cp IQTREE_results/covid_136_aligned.fasta.treefile covid_136_result_trees/covid_136_IQTree.nwk

# NeuralNJ: config/covid_136.yaml is the example config with instance_path "./examples/covid_136_aligned/".
# It needs ~890 GB, so it runs holding Phyloformer/.big_run.lock and only after the EB45
# Phyloformer boundary search has finished (that search's probe timer starts before flock).
cd NeuralNJ
python prepare_covid_phy.py 136         # -> examples/covid_136_aligned/{covid_136_aligned.phy,taxon_map.tsv}
flock ../Phyloformer/.big_run.lock python run_neuralnj_cpu.py finetune_rl_search.py \
    --config ./config/covid_136.yaml --infer_opt Argmax > NeuralNJ_results/covid_136_NeuralNJ.log 2>&1
cp output/Argmax_dim64_patch1/covid_136_aligned.tre NeuralNJ_results/
cp output/Argmax_dim64_patch1/covid_136_aligned.tre ../covid_136_result_trees/covid_136_NeuralNJ.nwk
cd ../Comparison
~/.conda/envs/Rphylo/bin/Rscript compare_simple_COVID_136.R
# -> figures_R/compare_simple_COVID_136.pdf, figures_R/quartets_common_COVID_136.csv
```

`compare_simple_COVID_136.R` is `compare_simple_COVID_300.R` with the covid_136 paths.

**Status (2026-10-29): the covid_136 NeuralNJ run has not happened.** It was queued to start after the
EB45 Phyloformer boundary search, but the queue script ran inside the Claude Code session and was
stopped when that session ended (before EB45 finished), so there is no covid_136 NeuralNJ log or tree.

**Update (2026-10-29 23:26):** raxmlpy passes the RAxMLpy README tests (`cd test; python test_raxmlpy.py`,
`test_pllpy_cmp_llh.py`, `test_pllpy_opt_br.py`). NeuralNJ cannot share the machine with the Phyloformer
boundary searches (`Phyloformer/run_pf_boundary_AS_CH.sh`): their probes reach 884-957 GB, leaving <= ~20 GB,
while NeuralNJ at full genome length needs > 58 GB even for 20 taxa (measured, stopped while still growing).
So covid_136 is queued with `NeuralNJ/run_covid_136_after_pf.sh`, started with
`setsid nohup NeuralNJ/run_covid_136_after_pf.sh > NeuralNJ/NeuralNJ_results/covid_136_queue.out 2>&1 &`:
it waits for the AS/CH searches to end, runs NeuralNJ holding `Phyloformer/.big_run.lock`, copies the tree to
`covid_136_result_trees/covid_136_NeuralNJ.nwk` and reruns `compare_simple_COVID_136.R`.

## 7. covid_10: first NeuralNJ tree on COVID (2026-10-30)

4 outgroups + 6 SARS-CoV-2 (`random.seed(42)`), 29,945 columns. Run while the AS45 Phyloformer search was
between probes, holding `Phyloformer/.big_run.lock` (AS45 probe 650 had finished; the next probe waited ~4 min).

```bash
cd COVID_data_fasta && python downsample_covid.py 10 && cd ..
CMAPLE/cmaple-1.1.0-Linux-intel/bin/cmaple -aln COVID_data_fasta/covid_10_aligned.fasta \
    -st DNA --prefix CMAPLE/CMAPLE_results/covid_10_aligned.fasta --seed 603
~/.conda/envs/iq_tree/bin/iqtree -s COVID_data_fasta/covid_10_aligned.fasta \
    -m GTR+G -T 4 --seed 603 --prefix IQTREE_results/covid_10_aligned.fasta
cd NeuralNJ
python prepare_covid_phy.py 10       # config/covid_10.yaml: instance_path "./examples/covid_10_aligned/"
setsid nohup ./run_covid_locked.sh 10 > NeuralNJ_results/covid_10_queue.out 2>&1 &
cd ../Comparison && ~/.conda/envs/Rphylo/bin/Rscript compare_simple_COVID_10.R
```

NeuralNJ writes the tree to `output/Argmax_dim64_patch1/covid_<N>_aligned.tre` (the trailing `/` in
`instance_path` makes the sub-folder name empty). The Argmax command has no `--branch_optimize`, so every
branch length is the placeholder `0.12345`; the comparison scripts drop them and treat the tree as topology only.

NeuralNJ: 4 min 7 s, peak 108 GB (`PEAK_KB 113122128`) - almost all of it the encoder's attention over
~30,000 sites, which does not depend on the number of taxa.

| covid_10 | IQ-TREE vs CMAPLE | IQ-TREE vs NeuralNJ |
|---|---|---|
| normalised RF | 0.000 | 0.692 |
| quartet distance (all 210 quartets) | 0.000 | 0.369 |
| IQ-TREE clades recovered | 100% | 50% |

Trees: `covid_10_result_trees/covid_10_{IQTree.nwk,CMAPLE.treefile,NeuralNJ.nwk}`;
figures: `Comparison/figures_R/compare_simple_COVID_10.pdf`.
