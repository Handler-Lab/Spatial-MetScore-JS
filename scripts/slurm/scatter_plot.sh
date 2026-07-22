#!/bin/bash
#SBATCH --job-name=plot_scatter
#SBATCH --output=logs/plot_scatter_%j.out
#SBATCH --error=logs/plot_scatter_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --partition=c128-m1024
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hhuan40@emory.edu
#
# Unified SLURM wrapper for the Spatial-MetScore Scatter Plotting.
# Note: Memory is reduced to 32G and time to 1 hour, which is more 
# than enough for reading a CSV and rendering ggplot2 figures.

set -euo pipefail

# Ensure logs directory exists
mkdir -p logs

# Navigate to working directory
cd /scratch/hhuan40/Spatial-MetScore-JS

# Load conda and activate the designated environment
source /opt/anaconda/etc/profile.d/conda.sh
conda activate /scratch/hhuan40/conda_envs/spatial_R_scratch

echo "=== Starting Scatter Plot Generation ==="
echo "Node: $(hostname)"
echo "Date: $(date)"

# Execute the R script (请确保这里是你保存画图代码的真实 R 脚本文件名)
Rscript scripts/scatter_plot.R

echo "=== Plotting Finished ==="
echo "Date: $(date)"