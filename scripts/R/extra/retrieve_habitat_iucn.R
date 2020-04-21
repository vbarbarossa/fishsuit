# RU CLUSTER --> could not install rredlist on cartesius

# retrieve habitat type
source('config.R');

libinv(c('sf','dplyr','rredlist'))

dir_habitat <- dir_(paste0('proc/iucn_habitat_list/'))

#my personal token
token <- 'd361026f05b472e57b0ffe1fa5c9a768aaf3d8391abbb464293e9efe2bbbf733'

# get unique species names
sp_names <- sf::read_sf('proc/species_ranges_merged.gpkg') %>% .$binomial

# retrieve the habitat type per species
for(i in sp_names[1:3]) {
  
  hab <- rl_habitats(as.character(i),key = token)$result$habitat
  saveRDS(hab,paste0(dir_habitat,i,'.rds'))
  
  Sys.sleep(2)
  
}