#!/bin/bash
# Re-check the AS45 (693/700) and EB45 (748/749) Phyloformer boundaries at low priority, so IQ-TREE on CH45_604 is not slowed.
# AS 693 peaked at 876 GB and EB 748 at 957 GB, so probes wait until other users hold < 120 GB (AS) / < 40 GB (EB).
cd "$(dirname "$0")/.."
pipeline/recheck_pf_boundary.sh AS 693 700 120
pipeline/recheck_pf_boundary.sh EB 748 750 40
