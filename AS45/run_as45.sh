#!/bin/bash
# Run every method on the AS45 "whole" supermatrix: 800 randomly selected loci (downsample_loci_800.py, seed 42).
# Phyloformer needs several hundred GB of RAM here, so it holds a lock and runs one dataset at a time.
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
ENVS=/mnt/home/jgao/.conda/envs
LOG=$ROOT/AS45/logs/runs
N=800
S=AS45_supermatrix_$N
mkdir -p $LOG
( cd $ROOT/IQTREE_results && $ENVS/iq_tree/bin/iqtree2 -s $S.fasta -m MFP -T AUTO --threads-max 28 -B 1000 -alrt 1000 --prefix $S ) > $LOG/iqtree.out 2>&1 &
( cd $ROOT/CMAPLE && ./cmaple-1.1.0-Linux-intel/bin/cmaple-aa -aln ../AS45_concat_$N/$S.fasta --prefix CMAPLE_results/$S -seed 697640 ) > $LOG/cmaple.out 2>&1 &
( cd $ROOT/MsaPhylo && $ENVS/msaphylo/bin/python MsaPhylo.py --i "./data/AS45_msas_$N/$S.fasta" --name AS45_$N --o "./AS45_embedding_trees_$N/" --l 2 ;
    $ENVS/msaphylo/bin/python run_as45_attention_trees_$N.py ) > $LOG/msaphylo.out 2>&1 &
( cd $ROOT/Phyloformer && flock $ROOT/Phyloformer/.big_run.lock bash -c "$ENVS/phylo/bin/python infer_alns.py -o data/testdata/as45_concat_matrices_$N models/pf.ckpt data/testdata/as45_msas_$N \
    && ./bin/bin_linux/fastme -i data/testdata/as45_concat_matrices_$N/$S.phy -o data/testdata/as45_concat_matrices_$N/$S.nwk --nni --spr" ) > $LOG/phyloformer.out 2>&1 &
wait
echo ALL DONE
