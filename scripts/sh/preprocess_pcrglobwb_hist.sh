#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH -p normal
#SBATCH --output=diag/preprocess_pcrglobwb_hist_year_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/pre/pcrglobwb_hist_year.R
