#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --array=1-17
#SBATCH --account=general
#SBATCH --partition=c64-m512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/logs/salmon_%A_%a.out
#SBATCH --error=/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/logs/salmon_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hhuan40@emory.edu

conda init bash > /dev/null 2>&1
source ~/.bashrc
conda activate st_env


WORK_DIR="/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/Reference/handler_2025"
CLEANED_DIR="${WORK_DIR}/cleaned"
INDEX_DIR="${WORK_DIR}/salmon_index_m33"
OUT_DIR="${WORK_DIR}/quantified"


mkdir -p ${OUT_DIR}

cd ${WORK_DIR}


SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SRR_Acc_List.txt)

echo "=============================================="
echo "Quantifying sample: $SRR"
echo "=============================================="


salmon quant -i ${INDEX_DIR} -l A \
    -1 ${CLEANED_DIR}/${SRR}_1_paired.fq.gz \
    -2 ${CLEANED_DIR}/${SRR}_2_paired.fq.gz \
    -p 4 -o ${OUT_DIR}/${SRR}_quant

echo "=============================================="
echo "Finished quantification for: $SRR"
echo "=============================================="