#!/bin/bash
# pipeline/stage3_from_probe.sh X N HI [THREADS] [WAIT_FOR_FAILED_N]
# Stage 3 at an N whose Phyloformer tree already exists from a boundary probe: reuse that tree (same seed-42 loci), run the
# rest, then collect, plot and refresh the workbook.  With WAIT_FOR_FAILED_N, first wait until run_stage3.sh reports
# that Phyloformer failed at that N, and do nothing if it succeeded instead.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=/mnt/home/jgao/.conda/envs/phylobench_ob45/bin/python
X=$1; N=$2; HI=$3; T=${4:-24}; WAIT=$5
x=$(echo $X | tr A-Z a-z)
LOG=$ROOT/${X}45/logs/stage3.log
PROBE=$ROOT/${X}45/logs/phyloformer_boundary/N$N/${X}45_supermatrix_$N.nwk
cd $ROOT

if [ -n "$WAIT" ]; then
  until grep -qE "Phyloformer (FAILED|tree done) at N=$WAIT" $LOG; do sleep 30; done
  grep -q "Phyloformer tree done at N=$WAIT" $LOG && exit 0
fi
[ -s $PROBE ] || { echo "$(date) no probe tree $PROBE" >> $LOG; exit 1; }

mkdir -p Phyloformer/data/testdata/${x}45_concat_matrices_$N
cp $PROBE Phyloformer/data/testdata/${x}45_concat_matrices_$N/
echo "$(date) start N=$N, Phyloformer tree reused from the boundary probe" >> $LOG
SKIP_PF=1 pipeline/run_subset.sh $X $N $T > ${X}45/logs/subset_${N}_all.out 2>&1 &
RUN=$!
sleep 60
COLS=$(sed -n 's/.*"supermatrix_length": \([0-9]*\).*/\1/p' ${X}45/logs/concat_summary_$N.json)
mkdir -p ${X}45/logs/phyloformer_boundary
printf 'largest success: %s loci (%s columns)\nsmallest failure above it: %s loci\n' $N $COLS $HI > ${X}45/logs/phyloformer_boundary/boundary.txt
pipeline/collect_result_trees.sh $X $N >> $LOG 2>&1

wait $RUN
pipeline/collect_result_trees.sh $X $N >> $LOG 2>&1
$PY pipeline/make_compare_R.py $X $N IQ,PF >> $LOG 2>&1
( source /etc/profile > /dev/null 2>&1; module load R-4.4.0; cd Comparison && Rscript compare_simple_${X}_${N}.R > figures_R/compare_simple_${X}_${N}.log 2>&1 )
flock $ROOT/pipeline/.summary.lock $PY pipeline/make_summary_xlsx.py >> $LOG 2>&1
echo "$(date) all done: Comparison/figures_R/compare_simple_${X}_${N}.pdf ($(grep identically Comparison/figures_R/compare_simple_${X}_${N}.log))" >> $LOG
