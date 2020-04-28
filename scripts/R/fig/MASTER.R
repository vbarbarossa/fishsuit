source('R/functions.R'); source('R/_dir_loc.R')

dir_mod <- dir_('fishsuit_completeRun_warming_4targets/')


climate_models <- c('gfdl','ipsl','hadgem','miroc','noresm')

warming_targets <- c('1.5','2.0','3.2','4.5')
warming_tab <- read.csv('GMT/thresholdYears_4targets.csv')

