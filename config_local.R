#----------------------------------------------------------------------------------
#>> model run name & location
model_name <- 'fishsuit'
model_folder <- ''

#----------------------------------------------------------------------------------
#>> climate models & scenarios settings
climate_models <- c('gfdl','ipsl','hadgem','miroc','noresm')
scenarios <- c('hist','rcp2p6','rcp4p5','rcp6p0','rcp8p5')
warming_targets <- as.character(format( c(1.5,2,3.2,4.5), nsmall = 1))

model_type <- 'warming_targets' # or 'year_targets'

# table for year thresholds
warming_tab <- read.csv('thresholdYears_4targets.csv')
colnames(warming_tab) <- c('model',warming_targets)

# threshold for minimum flow to consider when filtering the pcrglobwb maps
# in general a cross-sectional area in m2 = 2.77*Q^0.79 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003WR002082)
# then for a river of about 0.5 m2 section, flow has to be 0.1 m3/s
flow_filter_threshold <- 0.1

#----------------------------------------------------------------------------------
#>> DIRECTORIES
dir_master <- ''

dir_data <- paste0(dir_master,'data/')

dir_model <- ''

# source directory for PCR-GLOBWB output directories
# dir_src_pcrglobwb <- '/projects/0/milkun/'

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
source('scripts/R/fun/generic.R')


