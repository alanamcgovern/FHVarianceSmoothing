#!/bin/bash
#SBATCH --job-name=Compile
#SBATCH --partition short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1
#SBATCH -o Compile.out
#SBATCH -e Compile.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amcgov@uw.edu

module load R
export R_LIBS=~/R/4.3

Rscript Compile_Stan.R

