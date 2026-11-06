#!/bin/bash
# pipeline/collect_result_trees.sh X N [LABEL]
# Copy each method's tree for the <X>45 N-loci supermatrix into <X>45_<LABEL>_result_trees/ (LABEL defaults to N),
# with the same file names as OB45_*_result_trees.  Methods that failed or are still running are listed as MISSING;
# re-run once they finish.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
X=$1; N=$2; LABEL=${3:-$N}
x=$(echo $X | tr A-Z a-z)
S=${X}45_supermatrix_$N
OUT=$ROOT/${X}45_${LABEL}_result_trees
mkdir -p $OUT

copy() {  # copy <source> <destination>
  if [ -s "$1" ]; then cp "$1" "$2" && echo "  ok      $(basename "$2")"
  else echo "  MISSING $(basename "$2")  ($1)"; fi
}

echo "== $OUT"
copy $ROOT/IQTREE_results/$S.treefile                                 $OUT/${X}45_supermatrix_${LABEL}_IQTree.treefile
copy $ROOT/Phyloformer/data/testdata/${x}45_concat_matrices_$N/$S.nwk $OUT/${X}45_supermatrix_${LABEL}_Phyloformer.nwk
copy $ROOT/CMAPLE/CMAPLE_results/$S.treefile                          $OUT/${X}45_supermatrix_${LABEL}_CMAPLE.treefile
copy $ROOT/MsaPhylo/${X}45_embedding_trees_$N/${X}45_${N}_2.nwk       $OUT/${X}45_${LABEL}_2_embedding_tree.nwk
copy $ROOT/MsaPhylo/${X}45_attention_trees_$N/$S/${S}_0_4.nwk         $OUT/${X}45_supermatrix_${LABEL}_0_4_MSATransformer_AttentionTree.nwk
