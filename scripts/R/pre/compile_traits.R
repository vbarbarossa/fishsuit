# HABITAT TYPE FROM FISHBASE-----------------------------------------------------------------------------------------

source('config.R'); # always load as functions are loaded within this script

libinv(c('sf','dplyr','rfishbase'))
# set fishbase version
options(FISHBASE_VERSION="19.04")


df <- read_sf('proc/species_ranges_merged.gpkg')

# variables to retrieve:
#   initial range size (in km2)
#   habitat type
#   body length (in cm)
#   trophic group
#   commercial importance
#   climate zone

# RANGE SIZE ---------------------------------------------------------------------------------------------------------
df$area <- units::drop_units(st_area(st_transform(df,54009))/10**6)

# LOTIC-LENTIC HABITAT -----------------------------------------------------------------------------------------------

# use both IUCN and fishbase data to assign lotic-lentic habitat type
hiucn <- read.csv('proc/iucn_habitat_type.csv') %>%
  as_tibble %>%
  select(binomial,lotic,lentic)
# hiucn$lentic_only <- 0
# hiucn$lentic_only[hiucn$lentic == 1 & (hiucn$lotic == 0 | is.na(hiucn$lotic))] <- 1

hfishbase <- ecology(df$binomial) %>%
  select(binomial = Species, lotic = Stream, lentic = Lakes) %>%
  distinct(binomial,.keep_all = T)
hfishbase[hfishbase == -1] <- 1


hab <- bind_rows(
  # make first piece of table based on iucn data 
  data.frame(
    binomial = df$binomial[df$binomial %in% hiucn$binomial]
  ) %>% as_tibble() %>%
    left_join(hiucn)
  ,
  # and remaining species basedon fishbase
  data.frame(
    binomial = df$binomial[!df$binomial %in% hiucn$binomial]
  ) %>% as_tibble() %>%
    left_join(hfishbase)
) %>%
  arrange(binomial)
hab$lentic_only <- 0
hab$lentic_only[hab$lentic == 1 & (hab$lotic == 0 | is.na(hab$lotic))] <- 1

sum(hab$lentic_only) #1160

# FISHBASE TRAITS ----------------------------------------------------------------------------------------------------
# length, trophic, commercial importance, marine
fb1 <- species(df$binomial) %>%
  select(binomial = Species, length = Length, importance = Importance, marine = Saltwater) %>%
  distinct(binomial,.keep_all = T)
fb2 <- ecology(df$binomial) %>%
  select(binomial = Species, FoodTroph) %>%
  distinct(binomial,.keep_all = T)
  

fb2$FoodTroph[is.na(fb2$FoodTroph)] <- 0
fb2$foodtrophcat[fb2$FoodTroph > 0] <- 'Herbi.'
fb2$foodtrophcat[fb2$FoodTroph > 2.19 & fb2$FoodTroph <= 2.79] <- 'Omni.'
fb2$foodtrophcat[fb2$FoodTroph > 2.79] <- 'Carni.'
fb2$FoodTroph[fb2$FoodTroph == 0] <- NA

fb1$marine[fb1$marine == -1] <- 1

fb <- inner_join(fb1,fb2)

apply(fb,2,function(x) sum(is.na(x)))
# > nrow(fb)
# [1] 12934
# > apply(fb,2,function(x) sum(is.na(x)))
# binomial       length   importance    FoodTroph foodtrophcat 
# 0          753         9960         9346         9346 

# CLIMATE ZONES ------------------------------------------------------------------------------------------------------
cz <- read.csv('proc/climate_zones.csv')

# IUCN code ------------------------------------------------------------------------------------------------------
#iucn data for threat status**************************************************
token <- 'd361026f05b472e57b0ffe1fa5c9a768aaf3d8391abbb464293e9efe2bbbf733'
library(rredlist)
iucn_code <- foreach(ts = c("DD", "LC", "NT", "VU", "EN","CR", "EW", "EX", "LRlc", "LRnt", "LRcd"),.combine = 'rbind') %do%{
  t <- rl_sp_category(ts,key = token)$result %>% 
    as_tibble() %>%
    mutate(code = ts) %>%
    select(binomial = scientific_name,code)
  return(t)
} %>% arrange(binomial)

# MERGE TABLES --------------------------------------------------------------------------------------

tab <- df %>%
  as_tibble() %>%
  select(-geom) %>%
  left_join(hab) %>%
  left_join(fb) %>%
  left_join(cz) %>%
  left_join(iucn_code)

write.csv(tab,'proc/species_traits.csv',row.names = F)

