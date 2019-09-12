#sbatch --array=1-5
# one for each rcp*climate model
g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source('config.R'); # always load as functions are loaded within this script

libinv(c('raster','foreach'))

if(model_type == 'year_targets'){
  
  attr.tab <- expand.grid(climate_models,scenarios[scenarios!='hist'])
  
  clmod <- attr.tab[g,1]
  scen <- attr.tab[g,2]
  
  dir_pcrglobwb_out <- dir_(paste0(dir_model,clmod,'/pcrglobwb_processed/'))
  
  dir_src <- paste0(dir_src_pcrglobwb,'output_',clmod,'_',scen,'_mergetime/')
  
  # -------------
  
  # assign time span
  if(scen == 'hist'){
    years <- timespan_hist
    start_year <- 1951 #based on pcrglobwb output files
  }else{
    years <- timespan_scen
    start_year <- 2006 #based on pcrglobwb output files
  }
  
  # number of weeks in a year
  no.weeks <- 52
  
  # to make it more dynamic, later on can call an extra script 
  # where all the formulas for different metrics are store 
  # and compute only metrics requested by the user
  
  # metrics
  # mi: annual minimum weekly flow
  # zf: zero-flow weeks
  # cv: coefficient of variation
  # ff: flood frequncy
  varsQ <- c('Qmi','Qzf','Qcv','Qav')#'ff' once available
  varsT <- c('Tma','Tcv')
  
  calc_metrics <- function(x){
    
    area <- areas[x]
    
    dir_out <- dir_( # create areas directories inside clmod rirectory
      paste0(dir_pcrglobwb_out,'/',area,'/' )
    )
    
    Qfile <- list.files(paste0(dir_src,area,'/netcdf/'),pattern = 'discharge_weekAvg_output')
    Tfile <- list.files(paste0(dir_src,area,'/netcdf/'),pattern = 'waterTemp_weekAvg_output')
    
    # create list of brick-metrics to store the output of each year
    brickQ <- brick(raster(paste0(dir_src,area,'/netcdf/',Qfile),band=1))
    brickT <- brick(raster(paste0(dir_src,area,'/netcdf/',Tfile),band=1))
    res <- foreach(i = varsQ) %do% brickQ
    res <- c(res,foreach(i = varsT) %do% brickT)
    names(res) <- c(varsQ,varsT)
    
    # loop thorugh years and calculate the metrics
    for(yr in years){
      
      # need the index for the brick of each metric
      brickIndex <- yr - min(years) + 1
      
      # read one year as raster brick of 52 layers
      # set to zero cells with flow lower than 0.001
      rQ <- calc(
        brick(paste0(dir_src,area,'/netcdf/',Qfile))[[((yr-start_year)*no.weeks+1):((yr-start_year+1)*no.weeks)]],
        fun=function(x){x[x <= 0.001] <- 0; return(x)}
      )
      
      rT <- calc(
        brick(paste0(dir_src,area,'/netcdf/',Tfile))[[((yr-start_year)*no.weeks+1):((yr-start_year+1)*no.weeks)]],
        fun=function(x){x[x == 0 | x > (273.15+60)] <- NA; return(x)}
      )
      
      # calc metrics
      res[['Qav']][[brickIndex]] <- mean(rQ,na.rm = T)
      
      res[['Qmi']][[brickIndex]] <- min(rQ,na.rm = T)
      
      res[['Qzf']][[brickIndex]] <- sum(
        calc(rQ, fun=function(x){ x[x == 0] <- 1; x[x != 1] <- NA; return(x)} )
        ,na.rm=T)
      
      res[['Qcv']][[brickIndex]] <- calc(rQ,sd,na.rm = T)/calc(rQ,mean,na.rm=T)
      
      # calc metrics
      res[['Tma']][[brickIndex]] <- max(rT,na.rm = T)
      
      res[['Tcv']][[brickIndex]] <- calc(rT,sd,na.rm = T)/calc(rT,mean,na.rm=T)
      
      
    }
    
    # average over the 30 years
    res.av <- lapply(res,function(x) calc(x,mean,na.rm=T))
    
    # round up to 3 decimals (then automatically, values < 0.001 are set to zero)
    res.av <- lapply(res.av,function(x) round(x,3))
    res.av[['Qzf']] <- round(res.av[['Qzf']],0) # for Qzf
    
    # store results
    for(varname in c(varsQ,varsT)){
      writeRaster(res.av[[varname]],paste0(dir_out,varname,'_',scen,'.tif'),format="GTiff", overwrite=TRUE)
    }
    
    
  }
  
  parallel::mcmapply(calc_metrics,seq_along(areas),SIMPLIFY = F,mc.cores = ncores)
  
  
  # MERGE THE AREAS ---------------------------------------------------------------------
  
  dir_merged <- dir_(paste0(dir_pcrglobwb_out,'merged/'))
  
  metrics <- c(varsQ,varsT)
  
  Qavbin <- paste0(dir_merged,'Qavbin.tif')
  
  for(i in seq_along(metrics)){
    
    v <- list()
    for(a in seq_along(areas)){
      v[[a]] <- raster(paste0(dir_pcrglobwb_out,areas[a],'/',metrics[i],'_',scen,'.tif'))
    }
    
    names(v)[1:2] <- c('x', 'y')
    v$fun <- sum
    v$na.rm <- TRUE
    
    writeRaster(
      extend(do.call(mosaic, v),extent(raster(Qavbin)))*raster(Qavbin),
      paste0(dir_merged,metrics[i],'_',scen,'.tif'),
      format="GTiff", overwrite=TRUE
    )
    
  }
  
}

if(model_type == 'warming_targets'){
  
  library(foreach)
  attr.tab <- foreach(i = seq_along(warming_targets),.combine = 'rbind') %do% {
    data.frame(
      clmod = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[1]),
      scen = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[2]),
      year = warming_tab[,i+1],
      warmt = warming_targets[i])
  }
  attr.tab <- attr.tab[!is.na(attr.tab$year),]
  row.names(attr.tab) <- NULL
  
  clmod <- attr.tab[g,'clmod']
  scen <- attr.tab[g,'scen']
  year <- attr.tab[g,'year']
  warmt <- attr.tab[g,'warmt']
  
  dir_pcrglobwb_out <- dir_(paste0(dir_model,clmod,'/pcrglobwb_processed/'))
  
  dir_src <- paste0(dir_src_pcrglobwb,'output_',clmod,'_',scen,'_mergetime/')
  dir_src_hist <- paste0(dir_src_pcrglobwb,'output_',clmod,'_hist_mergetime/')
  
  # -------------
  
  years <- (year-14):(year+15)
  if(year == 2085) years <- (year-14):(year+14) # to adjust for the fact that time series ends at 2099 (and not at 2100)
  start_year <- 2006 #based on pcrglobwb output files
  start_year_hist <- 1951
  
  # number of weeks in a year
  no.weeks <- 52
  
  # to make it more dynamic, later on can call an extra script 
  # where all the formulas for different metrics are store 
  # and compute only metrics requested by the user
  
  # metrics
  # mi: annual minimum weekly flow
  # zf: zero-flow weeks
  # cv: coefficient of variation
  # ff: flood frequncy
  varsQ <- c('Qmi','Qzf','Qav')#'Qcv','ff' once available
  varsT <- c('Tma') #,'Tcv'
  
  calc_metrics <- function(x){
    
    area <- areas[x]
    
    dir_out <- dir_( # create areas directories inside clmod rirectory
      paste0(dir_pcrglobwb_out,'/',area,'/' )
    )
    
    Qfile <- list.files(paste0(dir_src,area,'/netcdf/'),pattern = 'discharge_weekAvg_output')
    Tfile <- list.files(paste0(dir_src,area,'/netcdf/'),pattern = 'waterTemp_weekAvg_output')
    
    # create list of brick-metrics to store the output of each year
    brickQ <- brick(raster(paste0(dir_src,area,'/netcdf/',Qfile),band=1))
    brickT <- brick(raster(paste0(dir_src,area,'/netcdf/',Tfile),band=1))
    res <- foreach(i = varsQ) %do% brickQ
    res <- c(res,foreach(i = varsT) %do% brickT)
    names(res) <- c(varsQ,varsT)
    
    # loop thorugh years and calculate the metrics
    for(yr in years){
      
      # need the index for the brick of each metric
      brickIndex <- yr - min(years) + 1
      
      if(yr < start_year){
        
        Qfileh <- list.files(paste0(dir_src_hist,area,'/netcdf/'),pattern = 'discharge_weekAvg_output')
        Tfileh <- list.files(paste0(dir_src_hist,area,'/netcdf/'),pattern = 'waterTemp_weekAvg_output')
        
        rQ <- calc(
          brick(paste0(dir_src_hist,area,'/netcdf/',Qfileh))[[((yr-start_year_hist)*no.weeks+1):((yr-start_year_hist+1)*no.weeks)]],
          fun=function(x){x[x <= 0.001] <- 0; return(x)}
        )
        
        rT <- calc(
          brick(paste0(dir_src_hist,area,'/netcdf/',Tfileh))[[((yr-start_year_hist)*no.weeks+1):((yr-start_year_hist+1)*no.weeks)]],
          fun=function(x){x[x == 0 | x > (273.15+60)] <- NA; return(x)}
        )
        
        
      }else{
        # read one year as raster brick of 52 layers
        # set to zero cells with flow lower than 0.001
        rQ <- calc(
          brick(paste0(dir_src,area,'/netcdf/',Qfile))[[((yr-start_year)*no.weeks+1):((yr-start_year+1)*no.weeks)]],
          fun=function(x){x[x <= 0.001] <- 0; return(x)}
        )
        
        rT <- calc(
          brick(paste0(dir_src,area,'/netcdf/',Tfile))[[((yr-start_year)*no.weeks+1):((yr-start_year+1)*no.weeks)]],
          fun=function(x){x[x == 0 | x > (273.15+60)] <- NA; return(x)}
        )
        
      }
      
      # calc metrics
      res[['Qav']][[brickIndex]] <- mean(rQ,na.rm = T)
      
      res[['Qmi']][[brickIndex]] <- min(rQ,na.rm = T)
      
      res[['Qzf']][[brickIndex]] <- sum(
        calc(rQ, fun=function(x){ x[x == 0] <- 1; x[x != 1] <- NA; return(x)} )
        ,na.rm=T)
      
      # res[['Qcv']][[brickIndex]] <- calc(rQ,sd,na.rm = T)/calc(rQ,mean,na.rm=T)
      
      # calc metrics
      res[['Tma']][[brickIndex]] <- max(rT,na.rm = T)
      
      # res[['Tcv']][[brickIndex]] <- calc(rT,sd,na.rm = T)/calc(rT,mean,na.rm=T)
      
      
    }
    
    # average over the 30 years
    res.av <- lapply(res,function(x) calc(x,mean,na.rm=T))
    
    # round up to 3 decimals (then automatically, values < 0.001 are set to zero)
    res.av <- lapply(res.av,function(x) round(x,3))
    res.av[['Qzf']] <- round(res.av[['Qzf']],0) # for Qzf
    
    # store results
    for(varname in c(varsQ,varsT)){
      writeRaster(res.av[[varname]],paste0(dir_out,varname,'_',scen,'_',warmt,'C_',year,'.tif'),format="GTiff", overwrite=TRUE)
    }
    
    
  }
  
  parallel::mcmapply(calc_metrics,seq_along(areas),SIMPLIFY = F,mc.cores = ncores)
  
  
  # MERGE THE AREAS ---------------------------------------------------------------------
  
  dir_merged <- dir_(paste0(dir_pcrglobwb_out,'merged/'))
  
  metrics <- c(varsQ,varsT)
  
  Qavbin <- paste0(dir_merged,'Qavbin.tif')
  
  for(i in seq_along(metrics)){
    
    v <- list()
    for(a in seq_along(areas)){
      v[[a]] <- raster(paste0(dir_pcrglobwb_out,areas[a],'/',metrics[i],'_',scen,'_',warmt,'C_',year,'.tif'))
    }
    
    names(v)[1:2] <- c('x', 'y')
    v$fun <- sum
    v$na.rm <- TRUE
    
    writeRaster(
      extend(do.call(mosaic, v),extent(raster(Qavbin)))*raster(Qavbin),
      paste0(dir_merged,metrics[i],'_',scen,'_',warmt,'C_',year,'.tif'),
      format="GTiff", overwrite=TRUE
    )
    
  }
  
  
}




