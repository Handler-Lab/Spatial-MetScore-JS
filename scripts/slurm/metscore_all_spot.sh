#!/bin/bash
#SBATCH --job-name=metscore
#SBATCH --output=logs/metscore_%j.out
#SBATCH --error=logs/metscore_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=c128-m1024
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hhuan40@emory.edu
#
# Unified SLURM wrapper for the MetScore pipeline.
# Note: Memory is set to 128G because ranking genes across the entire 
# combined spatial transcriptomics matrix is highly memory-intensive.

set -euo pipefail

# Ensure logs directory exists
mkdir -p logs

# Navigate to working directory
cd /scratch/hhuan40/Spatial-MetScore-JS

# Load conda and activate the designated environment
source /opt/anaconda/etc/profile.d/conda.sh
conda activate /group/jshandl-g00/spatial_R_scratch

echo "=== Starting MetScore Calculation ==="
echo "Node: $(hostname)"
echo "Date: $(date)"

# Execute the R script
Rscript scripts/metscore_all_spot.R

echo "=== MetScore Pipeline Finished ==="
echo "Date: $(date)"