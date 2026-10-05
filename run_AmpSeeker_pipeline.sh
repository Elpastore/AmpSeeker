#!/bin/bash
#SBATCH -A asadio
#SBATCH -J AmpSeeker
#SBATCH --output=logs/AmpSeeker.%J.out
#SBATCH --error=logs/AmpSeeker.%J.error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=asadio@mrc.gm
#SBATCH --time=28-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G

set -euo pipefail
mkdir -p logs

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ampseeker

cd /mnt/beegfs/projects/vector_malariagen_data/amplicon_data/AmpSeeker

# Clears any stale lock left by a previous run that died uncleanly (killed
# job, node reboot, timeout). Safe here because this job only runs once
# nothing else is queued/running against this directory — never run
# --unlock while another AmpSeeker job for this same directory is active
# (check `squeue -u asadio` first if submitting manually).
echo "${SLURM_CPUS_PER_TASK}"
snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda \
  --configfile config/eniyou_agam.yaml --unlock

snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda \
  --configfile config/eniyou_agam.yaml \
  --resources mem_mb=90000 \
  --rerun-incomplete --keep-going
