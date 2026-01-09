#!/bin/bash
#SBATCH --job-name=archive_pdff
#SBATCH --partition=interactive
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/archive_%j.out
#SBATCH --error=logs/archive_%j.err

# Get the current username
USERNAME=$(whoami)

BUILD="/mnt/parscratch/users/$USERNAME/data/iBEAt_Build/dixon-pdff"
ARCHIVE="/shared/abdominal_imaging/Archive/iBEAt_Build"
LOGIN_NODE="login1"

rsync -av --no-group --no-perms "$BUILD" "${LOGIN_NODE}:${ARCHIVE}"
