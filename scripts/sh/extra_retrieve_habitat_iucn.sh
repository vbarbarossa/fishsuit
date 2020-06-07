#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-20:00:00
#SBATCH -p normal
#SBATCH --output=diag/extra_retrieve_habitat_iucn.out
#SBATCH --mail-type=END
#SBATCH --mail-user=vbarbarossa@science.ru.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript scripts/R/xtr/retrieve_habitat_iucn.R
