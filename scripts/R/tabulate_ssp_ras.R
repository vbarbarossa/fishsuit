# ------------------------------------------------------------------------------------------------------------------------
#>> PREPROCESSING

source('config.R'); # always first

if(preprocess_iucn_data){
  
  libinv(c('rgdal','raster','sp'))
  
  preprocess_poly = FALSE
  preprocess_raster = TRUE
  preprocess_template = FALSE
  
  # create directories
  dir_points <- dir_(paste0(dir_proc,'ssp/ssp_points/'))
  dir_poly <- dir_(paste0(dir_proc,'ssp/ssp_poly/'))
  dir_ras <- dir_(paste0(dir_proc,'ssp/ssp_raster/'))
  
  # ------------------------------------------------------------------------------------------------------------------------
  #>> CREATE SINGLE SPECIES POLYGONS
  
  # save single species as single polygons
  if(preprocess_poly){
    
    # function to save single polygon species
    save_single_poly <- function(id) {
      
      saveRDS(sp_ranges[sp_ranges@data$id_no == id,],
              paste0(dir_poly,id,'.rds'))
      
    }
    
    # run the function in parallel for part 1 and 2 of the FW_FISH file
    #<<<<<<<< should switch to sf for faster loading of the shp data, but unfortunately could not install sf on cartesius
    for(part in 1:2){
      # read iucn shapefile
      sp_ranges <- readOGR(paste0(dir_data,'FW_FISH_20181113/FW_FISH_PART_',part,'.shp'),paste0('FW_FISH_PART_',part))
      
      # save single shapefiles to avoid overloading memory
      invisible(
        parallel::mcmapply(save_single_poly,as.character(sort(unique(sp_ranges@data$id_no))),mc.cores=ncores)
      )
      
      
    }
    
  }
  
  # read and save IDs
  ids <- sapply(strsplit(list.files(dir_poly,pattern = '.rds'),'\\.'),function(x) x[1])
  saveRDS(ids,paste0(dir_proc,'ssp/ids.rds'))
  
  # ------------------------------------------------------------------------------------------------------------------------
  #>> CREATE SINGLE SPECIES RASTER LAYERS AT 5 ARCMIN
  
  if(preprocess_raster){
    
    rasterize_gdal <- function(id,res=0.083333){
      
      dir_tmp <- dir_(paste0(dir_proc,'ssp/tmp',id,'/'))
      shape <- paste0(dir_tmp,id,'.shp')
      
      writeOGR(readRDS(paste0(dir_poly,id,'.rds')),
               shape,
               id,driver="ESRI Shapefile")
      
      ras <- paste0(dir_ras,id,'.tif')
      
      system(
        paste0(
          'gdal_rasterize -burn 1 -ot Byte -l ',id,' -te -180 -90 180 90 -tr ',res,' ',res,
          ' -co "COMPRESS=LZW" ',
          shape,' ',ras
        )
      )
      
      system(paste0('rm -r ',dir_tmp))
      
      
    }
    
    invisible(
      parallel::mcmapply(rasterize_gdal,ids,mc.cores = ncores)
    )
    
  }
  
  
  # create template xy table with all 5 arc min cells (based on the ldd of pcrglobwb)
  if(preprocess_template){
    
    t <- data.frame(rasterToPoints(raster(paste0(dir_data,'ldd.asc')),spatial = F)[,1:2])
    coordinates(t) <- c("x", "y")
    proj4string(t) <- proj4string(readRDS(paste0(dir_poly,'9.rds')))
    
    t <- as(t,'SpatialPointsDataFrame')
    t@data <- data.frame(occ = rep(NA,nrow(coordinates(t))))
    
    # save template point data
    saveRDS(t,paste0(dir_proc,'ssp/tab_template.rds'))
    
  }
  
  # ------------------------------------------------------------------------------------------------------------------------
  #>> TABULATE SINGLE SPECIES AT 5 ARC MIN
  
  tab_template <- readRDS(paste0(dir_proc,'ssp/tab_template.rds'))
  ids <- readRDS(paste0(dir_proc,'ssp/ids.rds'))
  
  # do only those not done yet
  ids_completed <- sapply(strsplit(list.files(dir_points,pattern = '.rds'),'\\.'),function(x) x[1])
  ids <- ids[!(ids %in% ids_completed)]
  
  tabulate_single_species <- function(id){
    
    tab <- tab_template
    
    # pick one species
    ssp <- raster(paste0(dir_ras,id,'.tif'))
    
    # overlay 
    tab@data$occ <- extract(ssp,tab)
    
    tab <- tab[tab@data$occ == 1,]
    
    # save as rds object (fastest and storage efficient)
    saveRDS(tab,paste0(dir_points,id,'.rds'))
    
  }
  
  # mc.cleanup = FALSE allows to keep on running new child processes even if some fail
  # this way can first compute smaller polygons from smaller areas on thin nodes
  # then use fat nodes to compute the more memory intensive
  parallel::mcmapply(tabulate_single_species,ids,mc.cleanup = FALSE,mc.cores = ncores)
  
}