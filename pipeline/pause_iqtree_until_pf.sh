#!/bin/bash
# pipeline/pause_iqtree_until_pf.sh PF_PID IQ_PID [IQ_PID ...]
# Pause IQ-TREE runs (no progress is lost) so a Phyloformer run gets the CPU, and resume them as soon as it exits.
PF=$1; shift
kill -STOP "$@"
trap 'kill -CONT "$@"' EXIT
while kill -0 $PF 2>/dev/null; do sleep 15; done
