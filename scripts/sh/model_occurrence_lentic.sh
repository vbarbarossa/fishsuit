#!/bin/bash
#SBATCH -N 1
#SBATCH -t 5:00:00
#SBATCH -p normal
#SBATCH --array 1-5
#SBATCH --output=diag/model_occurrence_lentic_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

Rscript scripts/R/extra/model_occurrence_lentic.R
