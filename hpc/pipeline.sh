#!/bin/bash
# Use the Bash shell to interpret the script

#SBATCH --partition=gpu         # GPU partition
#SBATCH --gres=gpu:1            # Request 1 GPU on that partition
#SBATCH --cpus-per-task=16      # Use 8 CPU cores for data loading
#SBATCH --mem=32G               # 32 GB of RAM
#SBATCH --time=72:00:00         # Up to 72 hours of training

# The cluster will send an email to this address if the job fails or ends
#SBATCH --mail-user=s.sourbron@sheffield.ac.uk
#SBATCH --mail-type=FAIL,END

# Assigns an internal “comment” (or name) to the job in the scheduler
#SBATCH --comment=dixon

# Assign a name to the job
#SBATCH --job-name=dixon

# Write logs to the logs folder
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Unsets the CPU binding policy.
# Some clusters automatically bind threads to cores; unsetting it can 
# prevent performance issues if your code manages threading itself 
# (e.g. OpenMP, NumPy, or PyTorch).
unset SLURM_CPU_BIND

# Ensures that all your environment variables from the submission 
# environment are passed into the job’s environment
export SLURM_EXPORT_ENV=ALL

# Loads the Anaconda module provided by the cluster.
# (On HPC systems, software is usually installed as “modules” to avoid version conflicts.)
module load Anaconda3/2024.02-1
module load Python/3.10.8-GCCcore-12.2.0 # essential to load latest GCC
module load CUDA/11.8.0 # must match with version in environment.yml

# Initialize Conda for this non-interactive shell
eval "$(conda shell.bash hook)"

# Activates your Conda environment named venv.
# (Older clusters use source activate; newer Conda versions use conda activate venv.)
# We assume that the conda environment has already been created
conda activate dixon

# Get the current username
USERNAME=$(whoami)

# Define path variables here
BUILD="/mnt/parscratch/users/$USERNAME/iBEAt_Build/dixon"

# Define a cache memory path explicitly so the default HOME dir is not used
# as this has limited space. This is to save model weights persistently
# and is also used for temporary storage of intermediate nifti files.
CACHE="/mnt/parscratch/users/$USERNAME/miblab-dl"

# Absolute path to the conda env you want to use (adjust if different)
ENV="/users/$USERNAME/.conda/envs/dixon"

# Prepend the environment's bin to PATH so subprocesses find CLI tools installed in the env
export PATH="$ENV/bin:$PATH"

# Ensure python can find installed site-packages (helpful if PYTHONPATH not set)
export PYTHONPATH="$ENV/lib/python3.10/site-packages:${PYTHONPATH:-}"

# Set CONDA_PREFIX so code that expects it can still work
export CONDA_PREFIX="$ENV"

export nnUNet_n_proc_DA=8 # Typically half the number of requested CPUs

# Diagnostics (will show in job output)
echo "---- Environment diagnostics (on compute node) ----"
echo "Using ENV: $ENV"
echo "Which python: $(which python || echo 'python not on PATH')"
echo "Explicit python: $ENV/bin/python"
echo "nnUNetv2_predict found at: $(which nnUNetv2_predict || echo 'NOT FOUND')"
echo "nnUNetv2_apply_postprocessing found at: $(which nnUNetv2_apply_postprocessing || echo 'NOT FOUND')"
echo "CONDA_PREFIX: $CONDA_PREFIX"
echo "LD_LIBRARY_PATH (start): ${LD_LIBRARY_PATH:-<empty>}"
echo "-----------------------------------------------"

# srun runs your program on the allocated compute resources managed by Slurm
srun "$ENV/bin/python" "$BUILD/iBEAt-dixon/src/stage_4_compute_fatwater.py" --build="$BUILD" --cache="$CACHE"