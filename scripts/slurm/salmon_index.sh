#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --account=general
#SBATCH --partition=c128-m1024
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/logs/salmon_index_%j.out
#SBATCH --error=/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/logs/salmon_index_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hhuan40@emory.edu


conda init bash > /dev/null 2>&1
source ~/.bashrc
conda activate st_env


cd /group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/Reference/handler_2025

echo "=============================================="
echo "Starting Salmon Indexing for M33 Reference..."
echo "=============================================="


salmon index -t gencode.vM33.transcripts.fa.gz -i salmon_index_m33 --gencode -p 8

echo "=============================================="
echo "Salmon Indexing Complete!"
echo "=============================================="