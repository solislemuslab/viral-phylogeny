#!/bin/bash
# CH45 stage-3 boundary to the nearest 10 below (600 ran; 800 failed above every success so far), probing only while
# other users hold < 60 GB of RAM; then run stage 3 at that N, reusing the probe's Phyloformer tree.
cd /mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
W=CH45/logs/phyloformer_boundary
PF_GATE_GB=60 /mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python pipeline/find_phyloformer_boundary.py CH 600 800 10 >> $W/search.out 2>&1
N=$(sed -n 's/largest success: \([0-9]*\) loci.*/\1/p' $W/boundary.txt)
HI=$(sed -n 's/smallest failure above it: \([0-9]*\) loci/\1/p' $W/boundary.txt)
[ "$N" = 600 ] && exit 0   # boundary unchanged: the existing CH45_600 stage-3 run stands
pipeline/stage3_from_probe.sh CH $N $HI 24
