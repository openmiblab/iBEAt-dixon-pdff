# iBEAt-dixon-pdff
PDFF calculation from DIXON data

On windows:

```bash
docker run --rm `
  -v "C:/Users/md1spsx/Documents/Data/iBEAt_Build:/data" `
  dixon_pdff `
  python src/stage_4_check_alignment.py --build /data
```

Deploy on HPC:

```bash
singularity build dixon_pdff.sif docker://dixon_pdff
```

Then run SLURM script

```bash
#!/bin/bash
#SBATCH --job-name=dixon_stage4
#SBATCH --output=stage4_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

module load singularity

# Get the current username
USERNAME=$(whoami)

BUILD="/mnt/parscratch/users/$USERNAME/data/iBEAt_Build"

singularity exec \
  -B "$BUILD:/data" \
  dixon_pdff.sif \
  python src/stage_4_check_alignment.py --build /data
```
