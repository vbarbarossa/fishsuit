#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH -p normal
#SBATCH --output=diag/preprocess_pcrglobwb_hist_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

Rscript scripts/R/preprocess_pcrglobwb_hist.R
