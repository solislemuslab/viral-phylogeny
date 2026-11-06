#!/bin/bash
# pipeline/auto_stage3.sh X [IQTREE_THREADS]
# Wait for <X>45/logs/phyloformer_boundary/boundary.txt (written by find_phyloformer_boundary.py), then run stage 3
# (IQ-TREE + Phyloformer on that many random loci), collect the trees, make the plots and refresh the summary workbook.
# Progress goes to <X>45/logs/auto_stage3.log.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=/mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python
X=$1; T=${2:-24}
B=$ROOT/${X}45/logs/phyloformer_boundary/boundary.txt
LOG=$ROOT/${X}45/logs/auto_stage3.log
cd $ROOT

echo "$(date) waiting for $B" >> $LOG
until [ -s $B ]; do sleep 300; done
N=$(sed -n 's/largest success: \([0-9]*\) loci.*/\1/p' $B)
echo "$(date) Phyloformer boundary N=$N; running stage 3" >> $LOG

pipeline/run_subset.sh $X $N $T > ${X}45/logs/subset_${N}_all.out 2>&1
echo "$(date) runs finished" >> $LOG

pipeline/collect_result_trees.sh $X $N >> $LOG 2>&1
METHODS=IQ
[ -s ${X}45_${N}_result_trees/${X}45_supermatrix_${N}_Phyloformer.nwk ] && METHODS=IQ,PF
$PY pipeline/make_compare_R.py $X $N $METHODS >> $LOG 2>&1
( source /etc/profile > /dev/null 2>&1; module load R-4.4.0; cd Comparison && Rscript compare_simple_${X}_${N}.R > figures_R/compare_simple_${X}_${N}.log 2>&1 )
echo "$(date) plots: Comparison/figures_R/compare_simple_${X}_${N}.pdf ($(grep identically Comparison/figures_R/compare_simple_${X}_${N}.log))" >> $LOG

# Only one copy should rewrite the workbook at a time
flock $ROOT/pipeline/.summary.lock $PY pipeline/make_summary_xlsx.py >> $LOG 2>&1
echo "$(date) done" >> $LOG
