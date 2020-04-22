#!/bin/bash
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH -p normal
#SBATCH --output=diag/preprocess_species2points.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/pre/species2points.R
