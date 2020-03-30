#!/bin/bash
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p normal
#SBATCH --array 1-5
#SBATCH --output=diag/model_occurrence_dispersal_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/extra/model_occurrence_dispersal.R
