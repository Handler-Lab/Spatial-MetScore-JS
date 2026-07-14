#!/bin/bash
#SBATCH --job-name=Spatial_MetScore
#SBATCH --output=logs/Spatial_MetScore%j.out
#SBATCH --error=logs/Spatial_MetScore%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=72:00:00
#SBATCH --partition=c128-m1024

# Send email notifications
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rprest2@emory.edu 

# Create logs directory if it doesn't exist
mkdir -p logs

echo "Starting Spatial MetScore Analysis..."

conda init bash > /dev/null 2>&1
source ~/.bashrc
conda activate spatial_R_scratch

# Run the Python script
Rscript /scratch/rprest2/Spatial-MetScore/scripts/copykat_pt1a.R

echo "Spatial MetScore Analysis complete."
