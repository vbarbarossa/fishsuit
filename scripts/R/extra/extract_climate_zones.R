
source('config.R'); # always load as functions are loaded within this script

library(raster); library(sf)

df <- read_sf('proc/species_ranges_merged.gpkg')

# read KG poly
kg <- raster('data/Koeppen-Geiger-Classification-Reclassfied_5min_moderesampling.tif')

# extract by polygon
list_val <- extract(kg,df)

# 
df <- df %>% as_tibble %>% dplyr::select(-geom)

df$climate_zone <- do.call('c',
        lapply(list_val,function(x){
          as.integer(names(sort(table(x),decreasing = T))[1])
        })
)

write.csv(df,'proc/climate_zones.csv',row.names = F)