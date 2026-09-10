#!/bin/bash
#SBATCH --job-name=Type2
#SBATCH --partition short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1-100
#SBATCH -o Type2.out
#SBATCH -e Type2.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amcgov@uw.edu

module load R
export R_LIBS=~/R/4.3

Rscript Type2_Models.R $1 "$SLURM_ARRAY_TASK_ID"

