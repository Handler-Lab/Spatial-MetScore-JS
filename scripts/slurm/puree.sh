#!/bin/bash
#SBATCH --job-name=Spatial_MetScore
#SBATCH --output=logs/Spatial_MetScore%j.out
#SBATCH --error=logs/Spatial_MetScore%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=c64-m512
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rprest2@emory.edu
#
# Unified SLURM wrapper for the PUREE pipeline (scripts/puree_pipeline.py).
# Resource envelope above is sized for the heaviest mode (perspot, all 54
# slides). Lighter modes (gmm, viz-*, pseudobulk) run fine within it, but you
# can shrink the request per-submission if the queue is busy, e.g.:
#   sbatch --mem=16G --cpus-per-task=4 --time=00:30:00 scripts/slurm/puree.sh gmm --k 3
#
# Usage:
#   sbatch scripts/slurm/puree.sh validate
#   sbatch scripts/slurm/puree.sh pseudobulk --output /path/to/puree_pseudobulk_54.csv
#   sbatch scripts/slurm/puree.sh perspot    --output /path/to/puree_perspot_allslides.csv
#   sbatch scripts/slurm/puree.sh viz-cohort
#   sbatch scripts/slurm/puree.sh viz-persample
#   sbatch scripts/slurm/puree.sh gmm --k 2
#   sbatch scripts/slurm/puree.sh gmm --k 3

set -euo pipefail
mkdir -p logs /group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/puree
cd /scratch/rprest2/Spatial-MetScore

CMD="${1:-}"; shift || true
case "$CMD" in
  validate|pseudobulk|perspot|gmm)
    PY=/scratch/rprest2/conda_puree/bin/python ;;   # needs sklearn / h5py / joblib
  viz-cohort|viz-persample)
    PY=/opt/anaconda/bin/python ;;                   # needs matplotlib
  *)
    echo "unknown or missing subcommand: '$CMD'" >&2
    echo "usage: sbatch puree.sh {validate|pseudobulk|perspot|viz-cohort|viz-persample|gmm} [args...]" >&2
    exit 1 ;;
esac

echo "=== $CMD ($PY) ==="
"$PY" scripts/puree_pipeline.py "$CMD" "$@"
echo DONE
