#!/bin/bash

#SBATCH --partition=long
#SBATCH --mem=36
#SBATCH --cpus-per-task=1
#SBATCH --job-name="get_eigenvecs"
#SBATCH -o get_eigenvecs_log.out
#SBATCH --mail-user=a.richard-bollans@kew.org
#SBATCH --mail-type=END,FAIL

set -euo # stop if fails
cd $SCRATCH/gentianales_trees/WCVP_12/Uphy
Rscript add_eigenvectors.R