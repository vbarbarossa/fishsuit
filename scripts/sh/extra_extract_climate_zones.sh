#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -p normal
#SBATCH --output=diag/extra_extract_climate_zones.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/extra/extract_climate_zones.R
