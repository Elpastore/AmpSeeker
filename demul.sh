#!/bin/bash
#SBATCH -A asadio
#SBATCH -J demultiplex
#SBATCH --error=demultiplex.%J.error
#SBATCH --output=demultiplex.%J.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=asadio@mrc.gm
#SBATCH --time=28-00:00:00
#SBATCH --ntasks=1



# Run the Python script
srun python ../build_multirun_metadata.py ../251211_M08382_0055_000000000-M783V  ../260618_M05061_0160_000000000-M78FW ../260622_M05061_0161_000000000-M78YC  ../260730_M05061_0163_000000000-M7837 ../260730_M08382_0056_000000000-M77NY --model config/example-metadata.tsv  --reads-dir ressources/cleaned_reads  --force-fastq -o config/metadata.tsv
