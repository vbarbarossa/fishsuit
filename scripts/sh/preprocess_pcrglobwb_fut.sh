#!/bin/bash
#SBATCH -N 1
#SBATCH -t 3:00:00
#SBATCH -p normal
#SBATCH --output=diag/preprocess_pcrglobwb_fut_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/pre/pcrglobwb_fut.R
