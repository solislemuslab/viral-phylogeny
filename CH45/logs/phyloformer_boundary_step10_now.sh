#!/bin/bash
# CH45 stage-3 boundary to the nearest 10 below, measured now under the current load (other users hold ~250 GB; each
# probe logs how much), so it is ready by 07:00 on 2026-11-06. 600 ran on 10-30 and 700 failed then, so the search runs between them (first probe 650).
# Then run stage 3 at that N, reusing the probe's Phyloformer tree.
cd /mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
W=CH45/logs/phyloformer_boundary
/mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python pipeline/find_phyloformer_boundary.py CH 600 700 10 >> $W/search.out 2>&1
N=$(sed -n 's/largest success: \([0-9]*\) loci.*/\1/p' $W/boundary.txt)
HI=$(sed -n 's/smallest failure above it: \([0-9]*\) loci/\1/p' $W/boundary.txt)
[ "$N" = 600 ] && exit 0   # boundary unchanged: the existing CH45_600 stage-3 run stands
pipeline/stage3_from_probe.sh CH $N $HI 24
