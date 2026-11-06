#!/bin/bash
# pipeline/run_subset.sh X N [THREADS]
# Build the <X>45 supermatrix of N random loci (pipeline/downsample_loci.py, seed 42) and run all five methods on it
# with the README commands: IQ-TREE, CMAPLE, MSA Transformer embedding + attention trees, Phyloformer (+ FastME).
# A method beyond its limit (MSA-T > 1023 columns, CMAPLE too divergent, Phyloformer out of memory) fails and says so in
# its log under <X>45/logs/subset_runs_<N>/.  THREADS is the IQ-TREE thread count (default 8).
# Set IQ_NICE=19 to give IQ-TREE low CPU priority, e.g. while a large Phyloformer run must finish first.
# Afterwards: pipeline/collect_result_trees.sh X N, pipeline/make_compare_R.py X N METHODS, Rscript the compare script.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
ENVS=/mnt/home/jgao/.conda/envs
X=$1; N=$2; T=${3:-8}
x=$(echo $X | tr A-Z a-z)
S=${X}45_supermatrix_$N
LOG=$ROOT/${X}45/logs/subset_runs_$N
mkdir -p $LOG

$ENVS/phylobench_ob45/bin/python $ROOT/pipeline/downsample_loci.py $X $N \
  && $ENVS/phylobench_ob45/bin/python $ROOT/pipeline/concatenate_supermatrix.py $X $N || exit 1
mkdir -p $ROOT/MsaPhylo/data/${X}45_msas_$N $ROOT/MsaPhylo/${X}45_embedding_trees_$N $ROOT/Phyloformer/data/testdata/${x}45_msas_$N
for d in $ROOT/IQTREE_results $ROOT/MsaPhylo/data/${X}45_msas_$N $ROOT/Phyloformer/data/testdata/${x}45_msas_$N; do
  cp $ROOT/${X}45_concat_$N/$S.fasta $d/
done

( cd $ROOT/IQTREE_results && ${IQ_NICE:+nice -n $IQ_NICE} $ENVS/iq_tree/bin/iqtree2 -s $S.fasta -m MFP -T $T -B 1000 -alrt 1000 --prefix $S ) > $LOG/iqtree.out 2>&1 &
( cd $ROOT/CMAPLE && ./cmaple-1.1.0-Linux-intel/bin/cmaple-aa -aln ../${X}45_concat_$N/$S.fasta --prefix CMAPLE_results/$S -seed 697640 ) > $LOG/cmaple.out 2>&1 &
( cd $ROOT/MsaPhylo && $ENVS/msaphylo/bin/python MsaPhylo.py --i "./data/${X}45_msas_$N/$S.fasta" --name ${X}45_$N --o "./${X}45_embedding_trees_$N/" --l 2 ;
    $ENVS/msaphylo/bin/python $ROOT/pipeline/run_attention_trees.py $X $N ) > $LOG/msaphylo.out 2>&1 &
# SKIP_PF=1: the Phyloformer tree already exists (e.g. from a boundary probe of the same N), so do not run it again
# Phyloformer needs ~5.5 GB of RAM per 1000 columns, so runs above 10k columns take turns through this lock
COLS=$(awk '/^>/ {if (n++) exit; next} {c += length($0)} END {print c}' $ROOT/IQTREE_results/$S.fasta)
LOCK=""; [ "$COLS" -gt 10000 ] && LOCK="flock $ROOT/Phyloformer/.big_run.lock"
[ -z "$SKIP_PF" ] && ( cd $ROOT/Phyloformer && $LOCK $ENVS/phylo/bin/python infer_alns.py -o data/testdata/${x}45_concat_matrices_$N models/pf.ckpt data/testdata/${x}45_msas_$N \
    && ./bin/bin_linux/fastme -i data/testdata/${x}45_concat_matrices_$N/$S.phy -o data/testdata/${x}45_concat_matrices_$N/$S.nwk --nni --spr ) > $LOG/phyloformer.out 2>&1 &
wait
echo ALL DONE
