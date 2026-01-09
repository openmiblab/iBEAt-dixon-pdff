#!/bin/bash
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=95:00:00
#SBATCH --mail-user=s.sourbron@sheffield.ac.uk
#SBATCH --mail-type=FAIL,END
#SBATCH --comment=align
#SBATCH --job-name=align
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

unset SLURM_CPU_BIND
export SLURM_EXPORT_ENV=ALL

module load Anaconda3/2024.02-1
module load Python/3.10.8-GCCcore-12.2.0 # essential to load latest GCC

# Get the current username
USERNAME=$(whoami)

# Define path variables here
ENV="/mnt/parscratch/users/$USERNAME/envs/pdff"
CODE="/mnt/parscratch/users/$USERNAME/scripts/iBEAt-dixon-pdff/src"
BUILD="/mnt/parscratch/users/$USERNAME/data/iBEAt_Build"
ARCHIVE="/shared/abdominal_imaging/Archive/iBEAt_Build"

# Initialize Conda for this non-interactive shell
# eval "$(conda shell.bash hook)"

# Activates your Conda environment named venv.
# (Older clusters use source activate; newer Conda versions use conda activate venv.)
# We assume that the conda environment has already been created
# conda activate "$ENV"

# srun runs your program on the allocated compute resources managed by Slurm
srun "$ENV/bin/python" "$CODE/pipeline.py" --build="$BUILD"