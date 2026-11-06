#!/bin/bash
# Exact CH45 stage-3 boundary between 600 (works) and 610 (fails), step 1: first probe 605, then bisection.
# If the answer is above 600, run stage 3 at that N, reusing the probe's Phyloformer tree.
cd /mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
W=CH45/logs/phyloformer_boundary
echo "$(date) exact search 600-610 step 1" >> $W/search.out
/mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python pipeline/find_phyloformer_boundary.py CH 600 610 1 >> $W/search.out 2>&1
N=$(sed -n 's/largest success: \([0-9]*\) loci.*/\1/p' $W/boundary.txt)
HI=$(sed -n 's/smallest failure above it: \([0-9]*\) loci/\1/p' $W/boundary.txt)
[ "$N" = 600 ] && exit 0
pipeline/stage3_from_probe.sh CH $N $HI 24
