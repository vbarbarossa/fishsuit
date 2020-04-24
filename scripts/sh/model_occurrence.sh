#!/bin/bash
#SBATCH -N 1
#SBATCH --array 1-5
#SBATCH -t 8:00:00
#SBATCH -p normal
#SBATCH --output=diag/model_occurrence_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/model_occurrence.R
