#!/bin/bash
#SBATCH --job-name=Combine
#SBATCH --partition short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1
#SBATCH -o Combine.out
#SBATCH -e Combine.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amcgov@uw.edu

module load R
export R_LIBS=~/R/4.3

Rscript Combine_Sims.R $1

