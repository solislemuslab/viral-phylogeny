#!/bin/bash
# Run every method on the shortest-loci OB45 supermatrices (15 and 20 loci).
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
ENVS=/mnt/home/jgao/.conda/envs
LOG=$ROOT/OB45/logs/short_runs
for N in 15 20; do
  S=OB45_supermatrix_short_$N
  ( cd $ROOT/IQTREE_results && $ENVS/iq_tree/bin/iqtree2 -s $S.fasta -m MFP -T AUTO -B 1000 -alrt 1000 --prefix $S ) > $LOG/iqtree_$N.out 2>&1 &
  ( cd $ROOT/CMAPLE && ./cmaple-1.1.0-Linux-intel/bin/cmaple-aa -aln ../OB45_concat_short_$N/$S.fasta --prefix CMAPLE_results/$S -seed 697640 ) > $LOG/cmaple_$N.out 2>&1 &
  ( cd $ROOT/MsaPhylo && $ENVS/msaphylo/bin/python MsaPhylo.py --i "./data/OB45_msas_short_$N/$S.fasta" --name OB45_short_$N --o "./OB45_embedding_trees_short_$N/" --l 2 \
      && $ENVS/msaphylo/bin/python run_ob45_attention_trees_short_$N.py ) > $LOG/msaphylo_$N.out 2>&1 &
  ( cd $ROOT/Phyloformer && $ENVS/phylo/bin/python infer_alns.py -o data/testdata/ob45_concat_matrices_short_$N models/pf.ckpt data/testdata/ob45_msas_short_$N \
      && ./bin/bin_linux/fastme -i data/testdata/ob45_concat_matrices_short_$N/$S.phy -o data/testdata/ob45_concat_matrices_short_$N/$S.nwk --nni --spr ) > $LOG/phyloformer_$N.out 2>&1 &
done
wait
echo ALL DONE
