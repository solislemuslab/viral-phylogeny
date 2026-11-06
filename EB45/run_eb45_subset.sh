#!/bin/bash
# run_eb45_subset.sh N [THREADS]: build the EB45 supermatrix of N random loci and run every method on it.
# Methods that exceed their limit (MSA Transformer > 1023 columns, CMAPLE too divergent, Phyloformer out of memory)
# fail and say so in their log, as on the full supermatrix.
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
ENVS=/mnt/home/jgao/.conda/envs
N=$1
T=${2:-8}
S=EB45_supermatrix_$N
LOG=$ROOT/EB45/logs/subset_runs_$N
mkdir -p $LOG

cd $ROOT/EB45 && $ENVS/phylobench_ob45/bin/python downsample_loci.py $N && $ENVS/phylobench_ob45/bin/python concatenate_supermatrix_n.py $N || exit 1
mkdir -p $ROOT/MsaPhylo/data/EB45_msas_$N $ROOT/MsaPhylo/EB45_embedding_trees_$N $ROOT/Phyloformer/data/testdata/eb45_msas_$N
for d in $ROOT/IQTREE_results $ROOT/MsaPhylo/data/EB45_msas_$N $ROOT/Phyloformer/data/testdata/eb45_msas_$N; do
  cp $ROOT/EB45_concat_$N/$S.fasta $d/
done

( cd $ROOT/IQTREE_results && $ENVS/iq_tree/bin/iqtree2 -s $S.fasta -m MFP -T $T -B 1000 -alrt 1000 --prefix $S ) > $LOG/iqtree.out 2>&1 &
( cd $ROOT/CMAPLE && ./cmaple-1.1.0-Linux-intel/bin/cmaple-aa -aln ../EB45_concat_$N/$S.fasta --prefix CMAPLE_results/$S -seed 697640 ) > $LOG/cmaple.out 2>&1 &
( cd $ROOT/MsaPhylo && $ENVS/msaphylo/bin/python MsaPhylo.py --i "./data/EB45_msas_$N/$S.fasta" --name EB45_$N --o "./EB45_embedding_trees_$N/" --l 2 ;
    $ENVS/msaphylo/bin/python run_eb45_attention_trees_n.py $N ) > $LOG/msaphylo.out 2>&1 &
# Phyloformer runs above 10k columns need hundreds of GB of RAM, so they take turns through this lock
COLS=$(awk '/^>/ {if (n++) exit; next} {c += length($0)} END {print c}' $ROOT/IQTREE_results/$S.fasta)
LOCK=""; [ "$COLS" -gt 10000 ] && LOCK="flock $ROOT/Phyloformer/.big_run.lock"
( cd $ROOT/Phyloformer && $LOCK $ENVS/phylo/bin/python infer_alns.py -o data/testdata/eb45_concat_matrices_$N models/pf.ckpt data/testdata/eb45_msas_$N \
    && ./bin/bin_linux/fastme -i data/testdata/eb45_concat_matrices_$N/$S.phy -o data/testdata/eb45_concat_matrices_$N/$S.nwk --nni --spr ) > $LOG/phyloformer.out 2>&1 &
wait
echo ALL DONE
