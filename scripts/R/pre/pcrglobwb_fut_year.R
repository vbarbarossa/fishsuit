#sbatch --array=1-5
# one for each rcp*climate model
g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source('config.R'); # always load as functions are loaded within this script

libinv(c('raster','foreach'))

attr.tab <- expand.grid(climate_models,scenarios[-1])

clmod <- as.character(attr.tab[g,1])
scen <- as.character(attr.tab[g,2])
year <- 2090
years <- 2080:2099
# warmt <- attr.tab[g,'warmt']

dir_pcrglobwb_out <- dir_(paste0(dir_model,'proc/additional_pcrglobwb_output/'))

dir_src <- paste0(dir_src_pcrglobwb,'output_',clmod,'_',scen,'_mergetime/')
dir_src_hist <- paste0(dir_src_pcrglobwb,'output_',clmod,'_hist_mergetime/')

# -------------

# years <- (year-14):(year+15)
# if(year == 2085) years <- (year-14):(year+14) # to adjust for the fact that time series ends at 2099 (and not at 2100)
start_year <- 2006 #based on pcrglobwb output files
start_year_hist <- 1951

# number of weeks in a year
no.weeks <- 52

# to make it more dynamic, later on can call an extra script 
# where all the formulas for different metrics are store 
# and compute only metrics requested by the user

# metrics
varsQ <- c('Qav')#'ff' once available
varsT <- c('Tav')

channel_section_area <- mask(raster(paste0(dir_data,'channel_section_area.tif')),raster(paste0(dir_data,'ldd.asc')) )
channel_section_area[channel_section_area[] < 0.1] <- NA

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
    # set to zero cells with flow lower than flow_filter_threshold
    rQ <- calc(
      brick(paste0(dir_src,area,'/netcdf/',Qfile))[[((yr-start_year)*no.weeks+1):((yr-start_year+1)*no.weeks)]],
      fun=function(x){x[x <= flow_filter_threshold] <- 0; return(x)}
    )
    
    rT <- calc(
      brick(paste0(dir_src,area,'/netcdf/',Tfile))[[((yr-start_year)*no.weeks+1):((yr-start_year+1)*no.weeks)]],
      fun=function(x){x[x == 0 | x > (273.15+60)] <- NA; return(x)}
    )
    
    
    # calc metrics
    res[['Qav']][[brickIndex]] <- mean(rQ,na.rm = T)
    res[['Tav']][[brickIndex]] <- mean(rT,na.rm = T)
    
    # res[['Qmi']][[brickIndex]] <- min(rQ,na.rm = T)
    
    # res[['Qma']][[brickIndex]] <- max(rQ,na.rm = T)
    
    # res[['Qve']][[brickIndex]] <- res[['Qma']][[brickIndex]]/mask(crop(channel_section_area,extent(rQ)),res[['Qma']][[brickIndex]])
    
    # res[['Qzf']][[brickIndex]] <- sum(
      # calc(rQ, fun=function(x){ x[x == 0] <- 1; x[x != 1] <- NA; return(x)} )
      # ,na.rm=T)
    
    # res[['Qcv']][[brickIndex]] <- calc(rQ,sd,na.rm = T)/calc(rQ,mean,na.rm=T)
    
    # calc metrics
    # res[['Tma']][[brickIndex]] <- max(rT,na.rm = T)
    
    # res[['Tmi']][[brickIndex]] <- min(rT,na.rm = T)
    
    # res[['Tcv']][[brickIndex]] <- calc(rT,sd,na.rm = T)/calc(rT,mean,na.rm=T)
    
    
  }
  
  # average over the 30 years
  res.av <- lapply(res,function(x) calc(x,mean,na.rm=T))
  
  # round up to 3 decimals (then automatically, values < flow_filter_threshold are set to zero)
  res.av <- lapply(res.av,function(x) round(x,3))
  # res.av[['Qzf']] <- round(res.av[['Qzf']],0) # for Qzf
  # res.av[['Tmi']] <- mask(res.av[['Tmi']],res.av[['Tmi']] > 273.2,maskvalue=0) # this is how mask works for binary filter layers
  
  # store results
  for(varname in c(varsQ,varsT)){
    writeRaster(res.av[[varname]],paste0(dir_out,varname,'_',clmod,'_',scen,'_',years[1],'-',years[length(years)],'.tif'),format="GTiff", overwrite=TRUE)
  }
  
  
}

parallel::mcmapply(calc_metrics,seq_along(areas),SIMPLIFY = F,mc.cores = ncores)


# MERGE THE AREAS ---------------------------------------------------------------------

dir_merged <- dir_(paste0(dir_pcrglobwb_out,'merged/'))

metrics <- c(varsQ,varsT)

# create binary layer based on long-term average annual flow
v <- list()
for(a in seq_along(areas)){
  v[[a]] <- raster(paste0(dir_pcrglobwb_out,areas[a],'/Qav_',clmod,'_',scen,'_',years[1],'-',years[length(years)],'.tif')) #<<<here should always be hist
}

names(v)[1:2] <- c('x', 'y')
v$fun <- sum
v$na.rm <- TRUE

Qavbin <- paste0(dir_merged,'Qavbin_',clmod,'_',scen,'_',years[1],'-',years[length(years)],'.tif')

l <- do.call(mosaic, v) >= flow_filter_threshold
NAvalue(l) <- 0 # set 0 to NA

writeRaster(l,Qavbin, format="GTiff", overwrite=TRUE)

Qavbin <- raster(Qavbin)
NAvalue(Qavbin) <- 0 # set 0 to NA

for(i in seq_along(metrics)){
  
  v <- list()
  for(a in seq_along(areas)){
    v[[a]] <- raster(paste0(dir_pcrglobwb_out,areas[a],'/',metrics[i],'_',clmod,'_',scen,'_',years[1],'-',years[length(years)],'.tif'))
  }
  
  names(v)[1:2] <- c('x', 'y')
  v$fun <- sum
  v$na.rm <- TRUE
  
  writeRaster(
    mask(extend(do.call(mosaic, v),extent(raster(Qavbin))),Qavbin),
    paste0(dir_merged,metrics[i],'_',clmod,'_',scen,'_',years[1],'-',years[length(years)],'.tif'),
    format="GTiff", overwrite=TRUE
  )
  
}

