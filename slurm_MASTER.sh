#!/bin/sh

cd /nfs/home2/valeriob/fishsuit_completeRun_warming_4targets/

job0=$(sbatch --parsable --array 1-5 scripts/sh/preprocess_pcrglobwb_hist.sh)

job1=$(sbatch --dependency=afterok:$job0 --parsable --array 1-42 scripts/sh/preprocess_pcrglobwb_fut.sh)

job2=$(sbatch --dependency=afterok:$job1 --parsable --array 1-30 scripts/sh/map_thresholds.sh)

job3=$(sbatch --dependency=afterok:$job2 --parsable --array 1-5 scripts/sh/model_occurrence.sh)
