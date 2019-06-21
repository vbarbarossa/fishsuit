#!/bin/bash
#SBATCH -N 1
#SBATCH -t 3:00:00
#SBATCH -p normal
#SBATCH --output=diag/preprocess_pcrglobwb_fut_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

Rscript scripts/R/preprocess_pcrglobwb_fut.R
