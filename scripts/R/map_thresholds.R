#sbatch --array=1-11

g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source('config.R'); # always load as functions are loaded within this script

libinv(c('raster','foreach'))

# array will be var1*var2 of attrs
attrs <- expand.grid(
  c('Qmi','Qma','Qzf','Qve','Tma','Tmi'), #'Qcv','Tcv'
  climate_models
)

metric <- as.character(attrs[g,1])
clmod <- as.character(attrs[g,2])

cat('\nclimate model = ',clmod,'\nmetric = ',metric,'\n')

yr <- c(timespan_hist[1],timespan_hist[2])

dir_merged <- paste0(dir_model,clmod,'/pcrglobwb_processed/merged/')
dir_niches <- dir_(paste0(dir_model,clmod,'/niches/'))

# read IUCN ids
iucn <- foreach(i = 1:2,.combine='rbind') %do% foreign::read.dbf(paste0(dir_data,'/FW_FISH_20181113/FW_FISH_PART_',i,'.dbf'))
iucn.id <- as.character(sort(unique(iucn$id_no)))

infiles <- list.files('proc/ssp/single_points',full.names = T)
ids <- lapply(infiles,function(x) strsplit(x[1],'/')[[1]][4] %>% strsplit(.,'\\.') %>% .[[1]] %>% .[1]) %>% do.call('c',.)

source('scripts/R/calc/map_variable2species.R')

map_variable2species(
  infile_lyr =  paste0(dir_merged,metric,'_hist.tif')#.tif
  ,infiles_pts = infiles #.rds # path to the point shapefiles
  ,ids = ids # vector of ID names
  ,outfile_name = paste0(dir_niches,metric) # name of output table (NO EXTENTION). The file is saved in csv and rds formats
  ,NC = ncores # number of cores to use
)

