#!/bin/bash
#SBATCH --job-name=Spatial_MetScore
#SBATCH --output=logs/Spatial_MetScore%j.out
#SBATCH --error=logs/Spatial_MetScore%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --partition=c64-m512
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rprest2@emory.edu
mkdir -p logs
cd /scratch/rprest2/Spatial-MetScore
/opt/anaconda/bin/python scripts/viz_copykat_grid.py
