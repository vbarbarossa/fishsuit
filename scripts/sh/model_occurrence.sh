#!/bin/bash
#SBATCH -N 1
#SBATCH -t 8:00:00
#SBATCH -p normal
#SBATCH --output=diag/model_occurrence_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

Rscript scripts/R/model_occurrence.R
