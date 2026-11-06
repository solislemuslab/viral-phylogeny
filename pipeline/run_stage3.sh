#!/bin/bash
# pipeline/run_stage3.sh X N HI [THREADS] [FALLBACK_N]
# Stage 3 at a chosen N: run all methods on N random loci, and as soon as Phyloformer succeeds record N as the
# Phyloformer boundary (HI = known failing N above it) and collect its tree.  If Phyloformer runs out of memory, stop this N
# and retry at FALLBACK_N (with HI = N) if one is given.  When IQ-TREE finishes: collect, plot, refresh the workbook.
# Progress goes to <X>45/logs/stage3.log.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=/mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python
X=$1; N=$2; HI=$3; T=${4:-24}; FB=$5
x=$(echo $X | tr A-Z a-z)
LOG=$ROOT/${X}45/logs/stage3.log
PFLOG=$ROOT/${X}45/logs/subset_runs_$N/phyloformer.out
TREE=$ROOT/Phyloformer/data/testdata/${x}45_concat_matrices_$N/${X}45_supermatrix_$N.nwk
cd $ROOT
rm -f $TREE $PFLOG

echo "$(date) start N=$N" >> $LOG
IQ_NICE=19 pipeline/run_subset.sh $X $N $T > ${X}45/logs/subset_${N}_all.out 2>&1 &
RUN=$!

until [ -s $TREE ] || grep -qE "RuntimeError|Traceback" $PFLOG 2>/dev/null; do sleep 30; done
if [ ! -s $TREE ]; then
  echo "$(date) Phyloformer FAILED at N=$N: $(grep -oE 'allocate [0-9]+ bytes' $PFLOG | tail -1)" >> $LOG
  pkill -u $USER -f "iqtree2 -s ${X}45_supermatrix_$N.fasta"
  pkill -u $USER -f "cmaple-aa -aln ../${X}45_concat_$N/"
  wait $RUN
  [ -n "$FB" ] && exec "$0" $X $FB $N $T
  exit 1
fi

COLS=$(sed -n 's/.*"supermatrix_length": \([0-9]*\).*/\1/p' ${X}45/logs/concat_summary_$N.json)
mkdir -p ${X}45/logs/phyloformer_boundary
printf 'largest success: %s loci (%s columns)\nsmallest failure above it: %s loci\n' $N $COLS $HI > ${X}45/logs/phyloformer_boundary/boundary.txt
pipeline/collect_result_trees.sh $X $N >> $LOG 2>&1
echo "$(date) Phyloformer tree done at N=$N ($COLS columns): ${X}45_${N}_result_trees/" >> $LOG

wait $RUN
pipeline/collect_result_trees.sh $X $N >> $LOG 2>&1
METHODS=IQ,PF
$PY pipeline/make_compare_R.py $X $N $METHODS >> $LOG 2>&1
( source /etc/profile > /dev/null 2>&1; module load R-4.4.0; cd Comparison && Rscript compare_simple_${X}_${N}.R > figures_R/compare_simple_${X}_${N}.log 2>&1 )
flock $ROOT/pipeline/.summary.lock $PY pipeline/make_summary_xlsx.py >> $LOG 2>&1
echo "$(date) all done: Comparison/figures_R/compare_simple_${X}_${N}.pdf ($(grep identically Comparison/figures_R/compare_simple_${X}_${N}.log))" >> $LOG
