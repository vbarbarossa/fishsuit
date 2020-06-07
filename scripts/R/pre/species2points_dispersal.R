
source('config.R'); # always load as functions are loaded within this script

# override, memory problems..
ncores <- 10

libinv(c('raster','foreach','dplyr','sf'))

# read feow_pcrglobwb watershed units
ws <- raster('data/BASINS_FEOW_PCRGLOBWB.tif')

# sample them with the template table
tab <- readRDS('proc/ssp/points_template.rds')
tab$ws <- extract(ws,tab)

dir_dispersal <- dir_('proc/ssp/single_points_dispersal/')

parallel::mcmapply(function(x){
  p <- readRDS(x)
  
  ws_ids <- tab[row.names(p),]$ws %>% unique
  
  saveRDS(tab[tab$ws %in% ws_ids,],paste0(dir_dispersal,tail(strsplit(x,'/')[[1]],n=1)))
  
},list.files('proc/ssp/single_points',full.names = T),mc.cores = ncores,SIMPLIFY = F)

