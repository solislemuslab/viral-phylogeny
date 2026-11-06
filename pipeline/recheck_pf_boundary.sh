#!/bin/bash
# pipeline/recheck_pf_boundary.sh X LO HI GATE_GB
# Re-check a Phyloformer boundary between LO (worked) and HI (failed) with step 1, starting each probe only when other
# users' processes hold less than GATE_GB of RAM (so the probe gets about as much memory as the original successes did).
# The existing boundary.txt is kept; the re-check result goes to boundary_recheck.txt and search_recheck.out.
cd "$(dirname "$0")/.."
X=$1; LO=$2; HI=$3; GATE=$4
W=${X}45/logs/phyloformer_boundary
cp $W/boundary.txt $W/boundary_before_recheck.txt
echo "$(date) recheck $LO-$HI step 1, gate: other users < $GATE GB" >> $W/search_recheck.out
PF_GATE_GB=$GATE /mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python pipeline/find_phyloformer_boundary.py $X $LO $HI 1 >> $W/search_recheck.out 2>&1
mv $W/boundary.txt $W/boundary_recheck.txt
cp $W/boundary_before_recheck.txt $W/boundary.txt
echo "$(date) done: $(head -1 $W/boundary_recheck.txt)" >> $W/search_recheck.out
