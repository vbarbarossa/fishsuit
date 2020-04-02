source('config.R'); # always load as functions are loaded within this script

libinv(c('sf','dplyr','foreach'))

# cat('\n\nCompiling datasets..\n')
# # AmazonFish dataset ------------------------------------------------------------
# amazonfish_tab <- read.csv2('data/amazonfish/CompleteDatabase.csv') %>% 
#   as_tibble %>%
#   filter(Occurrence.Status == 'valid') %>%
#   select(binomial = Referent.Species.Name, SBD_CD = SubDrainage.Code) %>%
#   mutate(binomial = gsub('.',' ',binomial,fixed=T)) %>%
#   distinct()
# 
# amazonfish <- read_sf('data/amazonfish/SubDrainageShapefile.shp') %>%
#   select(SBD_CD) %>%
#   right_join(amazonfish_tab) %>%
#   mutate(dataset = 'amazonfish') %>%
#   select(binomial,dataset) %>%
#   st_buffer(0)
# 
# amazonfish$area_sqkm = units::drop_units(st_area(st_transform(amazonfish,54009))/10**6)
# 
# 
# # IUCN dataset -----------------------------------------------------------------
# 
# iucn <- lapply(1:2,function(x) read_sf(paste0('data/FW_FISH_20200402/FW_FISH_PART',x,'.shp'))) %>%
#   do.call('rbind',.) %>%
#   filter(presence %in% 1:2) %>%
#   mutate(dataset = 'iucn') %>%
#   select(binomial,dataset)
# 
# iucn$area_sqkm = units::drop_units(st_area(st_transform(iucn,54009))/10**6)
# 
# 
# # PNAS paper dataset ------------------------------------------------------------------------------------------------
# cust <- read_sf('data/custom_ranges_poly.gpkg') %>%
#   filter(no_occ >= 10) %>%
#   mutate(dataset = 'custom') %>%
#   select(binomial = name, dataset) %>%
#   rename(geometry = geom)
# 
# cust$area_sqkm = units::drop_units(st_area(st_transform(cust,54009))/10**6)
# 
# 
# # OVERALL -----------------------------------------------------------------------------------------------------------
# sp <- do.call('rbind',list(amazonfish,iucn,cust)) %>%
#   select(binomial,area_sqkm,dataset)
# 
# 
# # synonyms check ----------------------------------------------------------------------------------------------------
# cat('Validating names in fishbase..\n')
# 
# # fishsuit names check
# library(rfishbase)
# options(FISHBASE_VERSION="19.04")
# 
# names <- sp$binomial %>% unique %>% as.character
# names_tab <- data.frame(binomial = names,fishbase_1 = NA, fishbase_2 = NA)
# for(i in seq_along(names)){
#   
#   n <- validate_names(names[i])
#   
#   if(length(n) > 0) names_tab[i,2] <- n[1]; names_tab[i,3] <- n[2] 
#   
# }
# 
# # some diagnostics
# 
# # number of species not in fishbase
# apply(names_tab,2,function(x) sum(!is.na(x)))
# # binomial fishbase_1 fishbase_2 
# # 13544      13148         30 
# 
# # number of unique names
# apply(names_tab,2,function(x) length(unique(x[!is.na(x)])))
# # binomial fishbase_1 fishbase_2 
# # 13544      12934         30 
# 
# # merge names
# sp_names <- left_join(sp,names_tab) %>%
#   select(binomial,fishbase_1,fishbase_2,area_sqkm, dataset)
# 
# # from which datasets are the NAs (names that could not be validated against fishbase)?
# sp_names %>% 
#   filter(is.na(fishbase_1)) %>%
#   select(binomial,dataset) %>%
#   as_tibble() %>%
#   select(-geometry) %>%
#   distinct() %>%
#   pull(dataset) %>% table()
# # amazonfish     custom       iucn 
# # 139          2        261 
# 
# 
# # total cumulative range area of species to exclude?
# (area_to_exclude <- sp_names %>% 
#     filter(is.na(fishbase_1)) %>%
#     pull(area_sqkm) %>% sum())
# # [1] 67520378 # 67M km2
# 
# # and in %?
# area_to_exclude/(sp %>% pull(area_sqkm) %>% sum) * 100
# # [1] 1.082322 # only 1%
# 
# # save the shapefile
# write_sf(sp_names,'proc/species_ranges_raw.gpkg')
# 
# # FILTERED ----------------------------------------------------------------------------------------------------------

cat('Compiling filtered dataset..\n')

sp_names <- read_sf('proc/species_ranges_raw.gpkg')

# filter out species without validated names against fishbase
sp_filtered <- sp_names %>%
  filter(!is.na(fishbase_1)) %>%
  select(binomial = fishbase_1,dataset)

# create an integer id
id_tab <- data.frame(binomial = sp_filtered$binomial %>% unique, id_no = 1:length(sp_filtered$binomial %>% unique))

# combine
sp_filtered <- left_join(sp_filtered,id_tab) %>%
  select(id_no,binomial,dataset)


# merge polygons toether ###
cat('Merging polygons..\n')

# define ids of species to merge
sp_to_merge = table(sp_filtered$binomial) %>% .[.>1] %>% names

to_merge <- sp_filtered %>% 
  filter(binomial %in% sp_to_merge)
to_merge <- split(to_merge,to_merge$id_no)

merged <- parallel::mcmapply(function(i) st_union(to_merge[[i]] %>% st_buffer(1/120)) %>% st_sf(),
                             seq_along(to_merge),mc.cores=10,SIMPLIFY=F) %>%
  do.call('rbind',.) %>%
  mutate(id_no = names(to_merge) %>% as.integer) %>%
  left_join(id_tab) %>%
  select(id_no,binomial)

merged$id_no <- names(to_merge)

sp_final <- rbind(
  merged,
  sp_filtered %>% filter(!binomial %in% sp_to_merge) %>% select(-dataset) %>% rename(geometry = geom)
)

cat('Saving..\n')

write_sf(sp_final,'proc/species_ranges_merged.gpkg')

# HABITAT TYPE FROM FISHBASE-----------------------------------------------------------------------------------------

# variables to retrieve:
#   initial range size (in km2)
#   body length (in cm)
#   climate zone
#   trophic group
#   habitat type
#   commercial importance

# RANGE SIZE


# retrieve list of freshwater fish species available from fishbase
# habitat <- ecology(custom_m$name) %>%
#   select(name = Species, Stream, Lake = Lakes) %>%
#   distinct()
# 
# habitat$OnlyLake <- 0
# habitat$OnlyLake[habitat$Lake == -1 & (habitat$Stream == 0 | is.na(habitat$Stream))] <- -1

