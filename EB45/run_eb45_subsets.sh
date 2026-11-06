#!/bin/bash
# Run every method on the shortest-loci EB45 supermatrices (26 and 46 loci).
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
ENVS=/mnt/home/jgao/.conda/envs
LOG=$ROOT/EB45/logs/subset_runs
mkdir -p $LOG
for N in 26 46; do
  S=EB45_supermatrix_$N
  ( cd $ROOT/IQTREE_results && $ENVS/iq_tree/bin/iqtree2 -s $S.fasta -m MFP -T AUTO --threads-max 8 -B 1000 -alrt 1000 --prefix $S ) > $LOG/iqtree_$N.out 2>&1 &
  ( cd $ROOT/CMAPLE && ./cmaple-1.1.0-Linux-intel/bin/cmaple-aa -aln ../EB45_concat_$N/$S.fasta --prefix CMAPLE_results/$S -seed 697640 ) > $LOG/cmaple_$N.out 2>&1 &
  ( cd $ROOT/MsaPhylo && mkdir -p ./EB45_embedding_trees_$N && $ENVS/msaphylo/bin/python MsaPhylo.py --i "./data/EB45_msas_$N/$S.fasta" --name EB45_$N --o "./EB45_embedding_trees_$N/" --l 2 ;
      $ENVS/msaphylo/bin/python run_eb45_attention_trees_$N.py ) > $LOG/msaphylo_$N.out 2>&1 &
  ( cd $ROOT/Phyloformer && $ENVS/phylo/bin/python infer_alns.py -o data/testdata/eb45_concat_matrices_$N models/pf.ckpt data/testdata/eb45_msas_$N \
      && ./bin/bin_linux/fastme -i data/testdata/eb45_concat_matrices_$N/$S.phy -o data/testdata/eb45_concat_matrices_$N/$S.nwk --nni --spr ) > $LOG/phyloformer_$N.out 2>&1 &
done
wait
echo ALL DONE
