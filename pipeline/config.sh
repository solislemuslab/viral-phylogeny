# Per-dataset settings shared by the pipeline scripts.  Source it, then call dataset_config <X>.
# OUTGROUP is the root tip of <X>45_reference_trees/<X>-astral.tre: every plot is rooted on it.
dataset_config() {
  case $1 in
    OB) OUTGROUP=MYROD ;;
    AG) OUTGROUP=METFS ;;
    AS) OUTGROUP=DOTSN ;;
    CH) OUTGROUP=DANRE ;;
    EB) OUTGROUP=CHLRE ;;
    *) echo "unknown dataset $1 (use OB, AG, AS, CH or EB)" >&2; return 1 ;;
  esac
}
