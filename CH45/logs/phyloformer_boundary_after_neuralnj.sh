#!/bin/bash
# Start the CH45 step-10 Phyloformer boundary search (600-700, first probe 650) only after NeuralNJ covid_300 has finished
# (user, 2026-11-05: never run it at the same time as NeuralNJ). "Finished" = a covid_300 tree or a covid_300 log ending
# (exit= / KILLED) newer than this script's start, and no NeuralNJ process left.  To start anyway: touch
# CH45/logs/phyloformer_boundary/GO
cd /mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia
W=CH45/logs/phyloformer_boundary
STAMP=$W/.waiting_since; touch $STAMP
echo "$(date) waiting for NeuralNJ covid_300 to finish (or for $W/GO)" >> $W/search.out
done300() {
  for f in NeuralNJ/output/Argmax_dim64_patch1/covid_300_aligned.tre NeuralNJ/NeuralNJ_results/covid_300_aligned.tre; do
    [ $f -nt $STAMP ] && return 0
  done
  L=NeuralNJ/NeuralNJ_results/covid_300_NeuralNJ.log
  [ $L -nt $STAMP ] && grep -qE "^exit=|KILLED" $L
}
neuralnj_running() { pgrep -u $USER -f "run_neuralnj_cpu.py|chain_200_then_300.sh|run_covid_locked.sh" > /dev/null; }
until [ -e $W/GO ] || { done300 && ! neuralnj_running; }; do sleep 60; done
echo "$(date) $([ -e $W/GO ] && echo 'GO file found' || echo 'NeuralNJ covid_300 finished'); starting the search" >> $W/search.out
rm -f $W/GO
exec CH45/logs/phyloformer_boundary_step10_now.sh
