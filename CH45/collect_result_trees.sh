#!/bin/bash
# Copy each method's CH45 result tree into CH45_<label>_result_trees/, using the same file names as OB45_*_result_trees.
# <label> is "short_<N>" for the N shortest loci and "<N>" for N random loci (as OB45_short_20 / OB45_7).
# Safe to re-run: trees that are not finished yet are reported as MISSING and picked up on the next run.
ROOT=/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
X=CH45

copy() {  # copy <source> <destination>
  if [ -s "$1" ]; then cp "$1" "$2" && echo "  ok      $(basename "$2")"
  else echo "  MISSING $(basename "$2")  ($1)"; fi
}

# label  run-name suffix (stage 1: all five methods, stage 2: no MSA-T, stage 3: IQ-TREE + Phyloformer only)
# (as in IQTREE_results/CH45_supermatrix<suffix>.*)
while read LABEL SUF; do
  S=${X}_supermatrix$SUF; PF=ch45_concat_matrices$SUF; TAG=$X$SUF
  OUT=$ROOT/${X}_${LABEL}_result_trees
  mkdir -p $OUT
  echo "== $OUT"
  copy $ROOT/IQTREE_results/$S.treefile                              $OUT/${X}_supermatrix_${LABEL}_IQTree.treefile
  copy $ROOT/Phyloformer/data/testdata/$PF/$S.nwk                    $OUT/${X}_supermatrix_${LABEL}_Phyloformer.nwk
  copy $ROOT/CMAPLE/CMAPLE_results/$S.treefile                       $OUT/${X}_supermatrix_${LABEL}_CMAPLE.treefile
  copy $ROOT/MsaPhylo/${X}_embedding_trees$SUF/${TAG}_2.nwk          $OUT/${X}_${LABEL}_2_embedding_tree.nwk
  copy $ROOT/MsaPhylo/${X}_attention_trees$SUF/$S/${S}_0_4.nwk       $OUT/${X}_supermatrix_${LABEL}_0_4_MSATransformer_AttentionTree.nwk
done <<EOF
5 _5
307 _307
EOF
