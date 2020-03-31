#!/bin/bash
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p broadwell_short
#SBATCH --output=diag/preprocess_species2points.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/preprocess_species2points.R
