#!/bin/bash
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH --array 1-8
#SBATCH -p normal
#SBATCH --output=diag/phyloreg_stoch_array_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/res/phyloreg_stochastic_array.R
