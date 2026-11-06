#!/bin/bash
# Run every method on the full AG45 supermatrix (286 loci).
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
ENVS=/mnt/home/jgao/.conda/envs
LOG=$ROOT/AG45/logs/runs
S=AG45_supermatrix
mkdir -p $LOG
( cd $ROOT/IQTREE_results && $ENVS/iq_tree/bin/iqtree2 -s $S.fasta -m MFP -T AUTO -B 1000 -alrt 1000 --prefix $S ) > $LOG/iqtree.out 2>&1 &
( cd $ROOT/CMAPLE && ./cmaple-1.1.0-Linux-intel/bin/cmaple-aa -aln ../AG45_concat/$S.fasta --prefix CMAPLE_results/$S -seed 697640 ) > $LOG/cmaple.out 2>&1 &
( cd $ROOT/MsaPhylo && $ENVS/msaphylo/bin/python MsaPhylo.py --i "./data/AG45_msas/$S.fasta" --name AG45 --o "./AG45_embedding_trees/" --l 2 ) > $LOG/msaphylo_embedding.out 2>&1 &
( cd $ROOT/MsaPhylo && $ENVS/msaphylo/bin/python run_ag45_attention_trees.py ) > $LOG/msaphylo_attention.out 2>&1 &
( cd $ROOT/Phyloformer && $ENVS/phylo/bin/python infer_alns.py -o data/testdata/ag45_concat_matrices models/pf.ckpt data/testdata/ag45_msas \
    && ./bin/bin_linux/fastme -i data/testdata/ag45_concat_matrices/$S.phy -o data/testdata/ag45_concat_matrices/$S.nwk --nni --spr ) > $LOG/phyloformer.out 2>&1 &
wait
echo ALL DONE
