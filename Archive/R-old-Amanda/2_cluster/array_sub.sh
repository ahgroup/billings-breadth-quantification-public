#!/bin/bash
#SBATCH --job-name=DimensionTestArray
#SBATCH --partition=batch        
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=als80846@uga.edu  
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=5gb   
#SBATCH --time=48:00:00   
#SBATCH --output=array_%x.%j.out
#SBATCH --error=array_%x.%j.err
#SBATCH --array=1-12
 
cd $SLURM_SUBMIT_DIR

module load R/4.0.0-foss-2019b   

export OMP_NUM_THREADS=6

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" input.lst)

Rscript array_test.R $file

