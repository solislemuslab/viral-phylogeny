#!/bin/bash
# Copy each method's AG45 result tree into AG45_<N>_result_trees/, using the same file names as OB45_*_result_trees.
# Safe to re-run: trees that are not finished yet are reported as MISSING and picked up on the next run.
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia

copy() {  # copy <source> <destination>
  if [ -s "$1" ]; then cp "$1" "$2" && echo "  ok      $(basename "$2")"
  else echo "  MISSING $(basename "$2")  ($1)"; fi
}

for N in 18 36 286; do
  if [ $N = 286 ]; then S=AG45_supermatrix; PF=ag45_concat_matrices; TAG=AG45
  else S=AG45_supermatrix_$N; PF=ag45_concat_matrices_$N; TAG=AG45_$N; fi
  OUT=$ROOT/AG45_${N}_result_trees
  mkdir -p $OUT
  echo "== $OUT"
  copy $ROOT/IQTREE_results/$S.treefile                              $OUT/AG45_supermatrix_${N}_IQTree.treefile
  copy $ROOT/Phyloformer/data/testdata/$PF/$S.nwk                    $OUT/AG45_supermatrix_${N}_Phyloformer.nwk
  copy $ROOT/CMAPLE/CMAPLE_results/$S.treefile                       $OUT/AG45_supermatrix_${N}_CMAPLE.treefile
  # MsaPhylo only fits <= 1023 columns, i.e. the 18-loci supermatrix
  if [ $N = 18 ]; then
    copy $ROOT/MsaPhylo/AG45_embedding_trees_$N/${TAG}_2.nwk         $OUT/AG45_${N}_2_embedding_tree.nwk
    copy $ROOT/MsaPhylo/AG45_attention_trees_$N/$S/${S}_0_4.nwk      $OUT/AG45_supermatrix_${N}_0_4_MSATransformer_AttentionTree.nwk
  fi
done
