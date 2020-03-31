
range2table <- function(
  in_shapefile_ranges #path to shapefile with species ranges
  ,out_shapefile_multipoints #name of file to be saved
  ,out_dir_shapefile_multipoints # name of output directory where the file will be saved
  ,id_col = 'id_no' #id column within the shapefile. If an id has multiple entries then all the rows corresponding to the id will be sampled
  ,ext = c(-180,180,-90,90) #xmin, xmax, ymin, ymax
  ,res = 1/12 #resolution in degrees
  ,mask_lyr = NULL #mask to filter out NAs, e.g. marine areas, to speed up calculations
  ,template_points = NULL #path for template table should be saved .rds
  ,sparse_table = NULL #path for sparse table should be saved .rds
  ,dir_single_ranges = NULL #directory where single ranges as point shapefiles should be saved
){
  library(valerioUtils)
  libinv(c('raster','sf','dplyr'))
  
  # generate template table from empty raster---------------------------------------
  
  cat('Creating template points shapefile table..\n')
  
  # create template raster
  ras <- raster(res = res, ext = extent(ext))
  ras[] <- 1
  if(!is.null(mask_lyr)) ras <- mask(ras,raster(mask_lyr)) #filter based on mask
  
  # convert to xy data frame
  xytab <- as.data.frame(coordinates(ras))[!is.na(as.data.frame(ras)),]
  
  # convert to sf
  tab <- st_as_sf(xytab, coords = c('x', 'y'), crs = crs(ras))
  
  # save template table shapefile
  if(!is.null(template_points)){
    cat('Saving template points shapefile..\n')
    saveRDS(tab,template_points)
  } 
  
  #---------------------------------------------------------------------------------
  
  ids = ranges %>% pull(get(id_col))
  if(is.character(in_shapefile_ranges)){
    ranges <- read_sf(in_shapefile_ranges)
  }else{
    ranges <- in_shapefile_ranges
  }
  
  cat('Computing sparse table..\n')
  lst <- st_contains(ranges,tab,sparse = T)
  
  cat('Merging and saving ranges as MULTIPOINT feature shapefile..\n')
  pts <- lapply(1:length(lst),
         function(n) tab[lst[[n]],] %>% mutate(id_no = ids[n])  %>% group_by(id_no) %>% summarize()
  ) %>% do.call('rbind',.) %>% group_by(id_no) %>% summarize()
  st_write(pts,paste0(dir_(out_dir_shapefile_multipoints),'/',out_shapefile_multipoints,'.gpkg'))
  
  if(!is.null(dir_single_ranges)){
    cat('Saving single ranges as POINT shapefiles..\n')
    invisible(lapply(1:nrow(pts),
                     function(n) saveRDS(pts[n,],paste0(dir_(dir_single_ranges),'/',pull(pts[n,],id_no),'.rds')) 
    ))
  }
  
  if(!is.null(sparse_table)){
    cat('Saving sparse table..\n')
    saveRDS(lst,sparse_table)
  }
  
  # 
  # OLD STUFF:
  # 
  # if(multiple_poly_per_id == F){
  #   cat('Reading ranges shapefile..\n')
  #   if(is.character(in_shapefile_ranges)){
  #     ranges <- read_sf(in_shapefile_ranges)
  #   }else{
  #     ranges <- in_shapefile_ranges
  #   }
  #   
  #   
  # }else{
  #   cat('Reading ranges shapefile and merging polygons by ids..\n')
  #   libinv(c('dplyr'))
  #   if(is.character(in_shapefile_ranges)){
  #     ranges <- read_sf(in_shapefile_ranges)
  #   }else{
  #     ranges <- in_shapefile_ranges
  #   }
  #   
  #   ranges <- ranges %>% 
  #     # st_set_precision(100) %>% 
  #     group_by(get(id_col)) %>% 
  #     summarise()
  #   
  # }
  # 
  # cat('Computing sparse table..\n')
  # lst <- st_contains(ranges,tab,sparse = T)
  # if(!is.null(sparse_table)){
  #   cat('Saving sparse table..\n')
  #   saveRDS(lst,sparse_table)
  # }
  # 
  # if(!is.null(dir_single_ranges)){
  #   cat('Saving single ranges as point shapefiles..\n')
  #   invisible(lapply(1:length(lst),
  #                    function(n) saveRDS(tab[lst[[n]],],paste0(dir_(dir_single_ranges),'/',ids[n],'.rds')) 
  #   ))
  # }
  # 
  
  
  # need to figure a way to save it as multipoints, not trivial
  # mp <- st_multipoint(cbind(rep(1,length(lst[[1]])),as.matrix(xytab[lst[[1]],])))
  
  
}