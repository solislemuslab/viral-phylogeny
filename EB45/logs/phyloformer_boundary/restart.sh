#!/bin/bash
# Re-run the EB45 Phyloformer boundary search between 642 (success) and 749 (failure) once the AS45_113
# Phyloformer run, which was not under the memory lock, has finished (it took ~130 GB during probes 668/695).
cd /mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
while ps -u $USER -o args= | grep -q "^/mnt/home/jgao/.conda/envs/phylo/bin/python infer_alns.py -o data/testdata/as45_concat_matrices_113"; do sleep 30; done
exec /mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python find_phyloformer_boundary.py EB 642 749
