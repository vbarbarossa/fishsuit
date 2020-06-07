source('config.R')

# copy R and sh folder into the new model folder
dir_model <- dir_(dir_model)

dir_diagnostics <- dir_(paste0(dir_model,'diag/'))

# copy scripts into new folder
system(paste0('cp -r scripts ',dir_model))

if(model_type == 'year_targets'){
  # create MASTER.sh
  cat(
    '#!/bin/sh',
    '\n\n',
    'cd ',dir_model,
    '\n\n',
    'job0=$(sbatch --parsable --array 1-',length(climate_models),' scripts/sh/preprocess_pcrglobwb_hist.sh)',
    '\n\n',
    'job1=$(sbatch --dependency=afterok:$job0 --parsable --array 1-',length(climate_models)*(length(scenarios)-1),' scripts/sh/preprocess_pcrglobwb_fut.sh)',
    '\n\n',
    'job2=$(sbatch --dependency=afterok:$job1 --parsable --array 1-',3*length(climate_models),' scripts/sh/map_thresholds.sh)', #3 variables*climate models
    '\n\n',
    'job3=$(sbatch --dependency=afterok:$job2 --parsable --array 1-',length(climate_models),' scripts/sh/model_occurrence.sh)',
    '\n',
    sep='',
    file = paste0(dir_model,'slurm_MASTER.sh'))
  
}

if(model_type == 'warming_targets'){
  
  # create MASTER.sh
  cat(
    '#!/bin/sh',
    '\n\n',
    'cd ',dir_model,
    '\n\n',
    'job0=$(sbatch --parsable --array 1-',length(climate_models),' scripts/sh/preprocess_pcrglobwb_hist.sh)',
    '\n\n',
    'job1=$(sbatch --dependency=afterok:$job0 --parsable --array 1-',sum(apply(warming_tab,1,function(x) {ncol(warming_tab)-1-sum(is.na(x))})),' scripts/sh/preprocess_pcrglobwb_fut.sh)',
    '\n\n',
    'job2=$(sbatch --dependency=afterok:$job1 --parsable --array 1-',3*length(climate_models),' scripts/sh/map_thresholds.sh)',
    '\n\n',
    'job3=$(sbatch --dependency=afterok:$job2 --parsable --array 1-',length(climate_models),' scripts/sh/model_occurrence.sh)',
    '\n',
    sep='',
    file = paste0(dir_model,'slurm_MASTER.sh'))
  
}

# copy config file
system(paste0('cp config.R ',dir_model))

# and year thresholds file
system(paste0('cp thresholdYears_4targets.csv ',dir_model))

# execute slurm_MASTER
system(paste0(
  'sh ',dir_model,'slurm_MASTER.sh'
))



