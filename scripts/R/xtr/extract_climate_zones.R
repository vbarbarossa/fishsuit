
source('config.R'); # always load as functions are loaded within this script

library(raster); library(sf); library(dplyr)

# override, memory problems..
ncores <- 10


# read KG poly
kg <- raster('data/Koeppen-Geiger-Classification-Reclassfied_5min_moderesampling.tif')

# sample them with the template table
tab <- readRDS('proc/ssp/points_template.rds')
tab$kg <- extract(kg,tab)

df <- parallel::mcmapply(function(x){
  p <- readRDS(x)
  
  kg_id <- tab[row.names(p),]$kg %>% table %>% sort(.,decreasing=T) %>% names %>% .[1]
  
  # in case the entire range does not have corresponding kg values
  # unlikely as the total NAs are about 0.4% for the template tab
  if(is.null(kg_id)) kg_id <- NA
  
  id <- tail(strsplit(x,'/')[[1]],n=1) %>% strsplit(.,'\\.') %>% .[[1]] %>% .[1]
  
  return(
    data.frame(
      id_no = id, climate_zone = kg_id
    )
  )
  
},list.files('proc/ssp/single_points',full.names = T),mc.cores = ncores,SIMPLIFY = F,USE.NAMES = F) %>%
  do.call('rbind',.)

write.csv(df,'proc/climate_zones.csv',row.names = F)

# OLD STUFF failed, too memory intensive for R..
# df <- read_sf('proc/species_ranges_merged.gpkg') %>%
#   filter(!st_is_empty(.))
# # extract by polygon
# list_val <- extract(kg,df)
# 
# # 
# df <- df %>% as_tibble %>% dplyr::select(-geom)
# 
# df$climate_zone <- do.call('c',
#         lapply(list_val,function(x){
#           as.integer(names(sort(table(x),decreasing = T))[1])
#         })
# )
