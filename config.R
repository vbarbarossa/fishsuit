#----------------------------------------------------------------------------------
#>> model run name & location
model_name <- 'fishsuit_completeRun_warming_4targets'
model_folder <- '/nfs/home2/valeriob/'

#----------------------------------------------------------------------------------
#>> climate models & scenarios settings
climate_models <- c('gfdl','ipsl','hadgem','miroc','noresm')
scenarios <- c('hist','rcp2p6','rcp4p5','rcp6p0','rcp8p5')
warming_targets <- as.character(format( c(1.5,2,3.2,4.5), nsmall = 1))

model_type <- 'warming_targets' # or 'year_targets'

# table for year thresholds
warming_tab <- read.csv('thresholdYears_4targets.csv')
colnames(warming_tab) <- c('model',warming_targets)

#----------------------------------------------------------------------------------
#>> DIRECTORIES
dir_master <- '/nfs/home2/valeriob/fishsuit/'

dir_data <- paste0(dir_master,'data/')

dir_model <- paste0(model_folder,model_name,'/')

# source directory for PCR-GLOBWB output directories
dir_src_pcrglobwb <- '/scratch-shared/jbosmans/'

#----------------------------------------------------------------------------------
#>> number of cores for parallelized scripts
ncores <- 22

#----------------------------------------------------------------------------------
#>> PCR-GLOBWB PRE-PROCESSING SETTINGS

# time spans for long-term averages (need to include all the years to average)
timespan_hist <- 1976:2005
timespan_scen <- 2036:2065 # for year_targets model

# areas from PCR-GLOBWB, remove any incomplete/defective
areas <- paste0('M',c(paste0(0,1:9),10:53))
areas <- areas[-29] # TEMPORARY remove the bugged M29

#----------------------------------------------------------------------------------
#>> IUCN DATA PRE-PROCESSING SETTINGS <<<<<<<<<< need to complete this part, keep it FALSE for now

preprocess_iucn_data <- FALSE


#----------------------------------------------------------------------------------
#>> MODEL OCCURRENCE
# variables
vars <- c('Qmi','Qzf','Tma')
thresholds <- c('2.5','97.5','97.5')

# minimum number of grid cells per species
min_no_grid_cells <- 10
filter_lentic_out <- FALSE

#----------------------------------------------------------------------------------
#>> CUSTOM FUNCTIONS
source('scripts/R/functions.R')


