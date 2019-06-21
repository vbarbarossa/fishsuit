#!/bin/bash
#SBATCH -N 1
#SBATCH -t 0:10:00
#SBATCH -p normal
#SBATCH --output=diag/map_thresholds_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

Rscript scripts/R/map_thresholds.R
