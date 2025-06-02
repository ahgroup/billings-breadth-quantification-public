#!/bin/bash
#SBATCH --job-name=BQ_main
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=6-00:00:00
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH -o logs/%x_%j.out

# Load necessary modules
module load R/4.4.2-gfbf-2024a

# Start targets pipeline
R -e "targets::tar_make()"
