#!/bin/bash
#SBATCH --job-name=copykat_all
#SBATCH --output=logs/copykat_%A_%a.out
#SBATCH --error=logs/copykat_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --partition=c128-m1024
#SBATCH --array=0-53%8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rprest2@emory.edu

# CopyKAT across all Visium samples, one array task per sample (max 8 concurrent).
cd /scratch/rprest2/Spatial-MetScore
mkdir -p logs
source /opt/anaconda/etc/profile.d/conda.sh
conda activate /group/jshandl-g00/spatial_R_scratch

MANIFEST=scripts/copykat_samples.txt
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MANIFEST")
echo "Array task $SLURM_ARRAY_TASK_ID -> sample $SAMPLE  (host $(hostname))"
Rscript scripts/copykat_sample.R "$SAMPLE"
echo "Task $SLURM_ARRAY_TASK_ID ($SAMPLE) finished with exit $?"
