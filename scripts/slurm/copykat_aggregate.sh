#!/bin/bash
#SBATCH --job-name=copykat_agg
#SBATCH --output=logs/copykat_agg_%j.out
#SBATCH --error=logs/copykat_agg_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --partition=c128-m1024
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rprest2@emory.edu

cd /scratch/rprest2/Spatial-MetScore
mkdir -p logs
source /opt/anaconda/etc/profile.d/conda.sh
conda activate /group/jshandl-g00/spatial_R_scratch
Rscript scripts/copykat_summarize.R
