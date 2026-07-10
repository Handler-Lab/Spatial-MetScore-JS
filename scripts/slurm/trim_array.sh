#!/bin/bash
#SBATCH --job-name=trim_array
#SBATCH --array=1-17
#SBATCH --account=general
#SBATCH --partition=c128-m1024
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/logs/trim_%A_%a.out
#SBATCH --error=/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/logs/trim_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hhuan40@emory.edu

conda init bash > /dev/null 2>&1
source ~/.bashrc
conda activate st_env

# Navigate to the raw data directory
cd /group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/Reference/handler_2025

# Fail-safe: Ensure the 'cleaned' directory exists in this folder before Trimmomatic runs
mkdir -p cleaned

# Extract the SRR ID for this specific array task
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SRR_Acc_List.txt)

echo "=============================================="
echo "Node ${SLURM_ARRAY_TASK_ID} is processing sample: $SRR"
echo "=============================================="

ADAPTER_PATH=$(find $CONDA_PREFIX/share -name "TruSeq3-PE.fa" | head -n 1)

# Execute trimming and output to the relative 'cleaned' folder
trimmomatic PE -threads 4 \
    ${SRR}_1.fastq ${SRR}_2.fastq \
    cleaned/${SRR}_1_paired.fq.gz cleaned/${SRR}_1_unpaired.fq.gz \
    cleaned/${SRR}_2_paired.fq.gz cleaned/${SRR}_2_unpaired.fq.gz \
    ILLUMINACLIP:${ADAPTER_PATH}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "=============================================="
echo "Node ${SLURM_ARRAY_TASK_ID} has finished processing sample: $SRR"
echo "=============================================="