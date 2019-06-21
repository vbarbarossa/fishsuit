#sbatch --array=1-11

slurm_arrayid<-Sys.getenv("SLURM_ARRAY_TASK_ID")
nodenr<-as.numeric(slurm_arrayid)

for(g in nodenr:nodenr){
  
  source('config.R'); # always load as functions are loaded within this script
  
  libinv(c('raster','foreach'))
  
  # array will be var1*var2 of attrs
  attrs <- expand.grid(
    c('Qmi','Qzf','Tma'), #'Qcv','Tcv'
    climate_models
  )
  
  metric <- as.character(attrs[g,1])
  clmod <- as.character(attrs[g,2])
  
  cat('\nclimate model = ',clmod,'\nmetric = ',metric,'\n')
  
  yr <- c(timespan_hist[1],timespan_hist[2])
  
  dir_merged <- paste0(dir_model,clmod,'/pcrglobwb_processed/merged/')
  dir_niches <- dir_(paste0(dir_model,clmod,'/niches/'))
  
  # set temporary directory for calculations
  dir_tmp <- dir_(paste0('tmp_nichesCalc_',metric,clmod,'_/'))
  rasterOptions(tmpdir = dir_tmp)
  
  # read IUCN ids
  iucn <- foreach(i = 1:2,.combine='rbind') %do% foreign::read.dbf(paste0(dir_data,'/FW_FISH_20181113/FW_FISH_PART_',i,'.dbf'))
  iucn.id <- as.character(sort(unique(iucn$id_no)))
  
  # load pcrglobwb discharge output
  Q <- raster(paste0(dir_merged,metric,'_hist.tif'))
  
  # function to retrieve the range quantiles
  range_custom <- function(id){
    
    pts <- readRDS(paste0(dir_master,'ssp/ssp_points/',id,'.rds'))
    
    # extract values
    val <- extract(Q,pts)
    val <- val[!is.na(val)]
    
    # retrieve quantiles of the values distribution
    q <- quantile(val,seq(0,1,0.005))
    
    # return the data as data.frame
    return(
      cbind(data.frame(IUCN_ID = id, no.grids = length(val), mean = mean(val), sd = sd(val)),
            t(as.data.frame(q)))
    )
    
  }
  
  # apply the function in parallel to the iucn ids
  t_range <- do.call(
    'rbind',
    parallel::mcmapply(range_custom,iucn.id,SIMPLIFY = F,mc.cores = ncores)
  )
  row.names(t_range) <- NULL
  
  # save results
  write.csv(t_range,paste0(dir_niches,metric,'.csv'), row.names = FALSE)
  saveRDS(t_range,paste0(dir_niches,metric,'.rds'))
  
  # remove temporary directory
  system(
    paste0('rm -r ',dir_tmp)
  )
  
}
