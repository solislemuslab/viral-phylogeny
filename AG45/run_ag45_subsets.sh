#!/bin/bash
# Run every method on the shortest-loci AG45 supermatrices (18 and 36 loci).
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
ENVS=/mnt/home/jgao/.conda/envs
LOG=$ROOT/AG45/logs/subset_runs
mkdir -p $LOG
for N in 18 36; do
  S=AG45_supermatrix_$N
  ( cd $ROOT/IQTREE_results && $ENVS/iq_tree/bin/iqtree2 -s $S.fasta -m MFP -T AUTO -B 1000 -alrt 1000 --prefix $S ) > $LOG/iqtree_$N.out 2>&1 &
  ( cd $ROOT/CMAPLE && ./cmaple-1.1.0-Linux-intel/bin/cmaple-aa -aln ../AG45_concat_$N/$S.fasta --prefix CMAPLE_results/$S -seed 697640 ) > $LOG/cmaple_$N.out 2>&1 &
  ( cd $ROOT/MsaPhylo && $ENVS/msaphylo/bin/python MsaPhylo.py --i "./data/AG45_msas_$N/$S.fasta" --name AG45_$N --o "./AG45_embedding_trees_$N/" --l 2 ;
      $ENVS/msaphylo/bin/python run_ag45_attention_trees_$N.py ) > $LOG/msaphylo_$N.out 2>&1 &
  ( cd $ROOT/Phyloformer && $ENVS/phylo/bin/python infer_alns.py -o data/testdata/ag45_concat_matrices_$N models/pf.ckpt data/testdata/ag45_msas_$N \
      && ./bin/bin_linux/fastme -i data/testdata/ag45_concat_matrices_$N/$S.phy -o data/testdata/ag45_concat_matrices_$N/$S.nwk --nni --spr ) > $LOG/phyloformer_$N.out 2>&1 &
done
wait
echo ALL DONE
