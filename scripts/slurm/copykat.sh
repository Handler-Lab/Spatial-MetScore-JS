#!/bin/bash
#SBATCH --job-name=copykat
#SBATCH --output=logs/copykat_%j.out
#SBATCH --error=logs/copykat_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --partition=c128-m1024
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rprest2@emory.edu
#
# Unified SLURM wrapper for the CopyKAT pipeline (scripts/copykat_pipeline.R).
# Resource envelope above is sized for a single-sample run (mode "sample").
# "combined" needs far more memory/time/fewer cores and should be submitted
# with explicit overrides, e.g.:
#   sbatch --mem=500G --cpus-per-task=4 --time=72:00:00 scripts/slurm/copykat.sh combined
# "summarize" is cheap and can shrink the request, e.g.:
#   sbatch --mem=8G --cpus-per-task=1 --time=00:20:00 scripts/slurm/copykat.sh summarize
#
# Usage:
#   sbatch scripts/slurm/copykat.sh sample <SAMPLE_ID> [--outdir DIR] [--distance TYPE] [--ncores N]
#   sbatch --array=0-53%8 scripts/slurm/copykat.sh sample   # (see note below on array mode)
#   sbatch --mem=500G --time=72:00:00 scripts/slurm/copykat.sh combined [--outdir DIR] [--ncores N]
#   sbatch scripts/slurm/copykat.sh summarize [--outdir DIR] [--out CSV]
#
# Array-job note: when submitted with --array, SAMPLE is NOT read from argv --
# it is looked up in scripts/copykat_samples.txt by $SLURM_ARRAY_TASK_ID, so
# the per-sample cohort sweep is:
#   sbatch --array=0-53%8 scripts/slurm/copykat.sh sample

set -euo pipefail
mkdir -p logs
cd /scratch/rprest2/Spatial-MetScore
source /opt/anaconda/etc/profile.d/conda.sh
conda activate /group/jshandl-g00/spatial_R_scratch

CMD="${1:-}"; shift || true

if [[ "$CMD" == "sample" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  MANIFEST=scripts/copykat_samples.txt
  SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MANIFEST")
  echo "Array task $SLURM_ARRAY_TASK_ID -> sample $SAMPLE  (host $(hostname))"
  Rscript scripts/copykat_pipeline.R sample "$SAMPLE" "$@"
  echo "Task $SLURM_ARRAY_TASK_ID ($SAMPLE) finished with exit $?"
  exit 0
fi

case "$CMD" in
  sample|combined|summarize)
    echo "=== $CMD ==="
    Rscript scripts/copykat_pipeline.R "$CMD" "$@"
    echo DONE ;;
  *)
    echo "unknown or missing subcommand: '$CMD'" >&2
    echo "usage: sbatch copykat.sh {sample <SAMPLE_ID>|combined|summarize} [args...]" >&2
    exit 1 ;;
esac
