# HABITAT TYPE FROM FISHBASE-----------------------------------------------------------------------------------------

source('config.R'); # always load as functions are loaded within this script

libinv(c('sf','dplyr','rfishbase'))
# set fishbase version
options(FISHBASE_VERSION="19.04")


df <- read_sf('proc/species_ranges_merged.gpkg')

# variables to retrieve:
#   initial range size (in km2)
#   body length (in cm)
#   climate zone
#   trophic group
#   habitat type
#   commercial importance

# RANGE SIZE ---------------------------------------------------------------------------------------------------------
df$area <- units::drop_units(st_area(st_transform(df,54009))/10**6)

# LOTIC-LENTIC -------------------------------------------------------------------------------------------------------

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


tab <- bind_rows(
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
tab$lentic_only <- 0
tab$lentic_only[tab$lentic == 1 & (tab$lotic == 0 | is.na(tab$lotic))] <- 1

# FISHBASE TRAITS ----------------------------------------------------------------------------------------------------


# CLIMATE ZONES ------------------------------------------------------------------------------------------------------
library(raster)
# read KG poly
kg <- raster('data/Koeppen-Geiger-Classification-Reclassfied_5min_moderesampling.tif')

extract(kg,df[1:5,]) %>% 
# extract by polygon


