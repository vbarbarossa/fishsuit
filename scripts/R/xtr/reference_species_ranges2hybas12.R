source('config.R')

# borrowed from connectfish

# determine centroids with sf
# if(file.exists('proc/hybas12_points_nolakes.gpkg')){
#   points <- read_sf('proc/hybas12_points_nolakes.gpkg')
# }else{
#   points <- foreach(cont = c('af','ar','as','au','eu','gr','na','sa','si'),.combine='rbind') %do% {
#     poly <- read_sf(paste0(dir_hybas12,'/hybas_',cont,'_lev12_v1c.shp'))
#     return(st_centroid(poly))
#   }
#   write_sf(points,'proc/hybas12_points_nolakes.gpkg',driver='GPKG')
# }


library(sf); library(dplyr)

points <- read_sf('data/hybas12_points_nolakes.gpkg')

# load species shapefile
sp <- read_sf('proc/species_ranges_raw.gpkg')

# reference to hydrobasins level 12
lst <- st_contains(sp,points,sparse = T)
# lst is a sparse matrix where each entry is a row of sp and contains a list of hybas12 points falling within that species polygos

# make database where each entry is a hybas ID and 
# loop through the species
# for each species, create a table with hybasID and species id_no

# should update this part using dplyr as in extract_customRanges2hybas12.R
tab <- lapply(seq_along(lst),function(i){
  hb <- points$HYBAS_ID[lst[[i]]]
  if(length(hb) > 0){
    return(
      data.frame(HYBAS_ID = hb,
                 binomial = sp$binomial[i])
    )
  }
}
) %>% do.call('rbind',.) %>% distinct()

saveRDS(tab,'proc/species_ranges_raw_on_hybas12.rds')
write.csv(tab,'proc/species_ranges_raw_on_hybas12.csv',row.names = F)
