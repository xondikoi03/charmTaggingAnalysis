#!/bin/bash
#SBATCH --job-name=WcAnalysis
#SBATCH --output=logs/WcAnalysis_%j.out
#SBATCH --error=logs/WcAnalysis_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=standard

# ===============================
# Environment setup
# ===============================

echo "Job started on $(hostname)"
echo "SLURM job ID: $SLURM_JOB_ID"

# Load modules if your cluster uses them
# module load python/3.10
# module load gcc

# Activate your environment
# ---- EDIT THIS ----
source ~/miniconda3/etc/profile.d/conda.sh
conda activate coffea-env

# Make sure matplotlib does not try to open X windows
export MPLBACKEND=Agg

# Optional: limit thread over-subscription
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ===============================
# Run analysis
# ===============================

python WcAnalysis.py

echo "Job finished at $(date)"
