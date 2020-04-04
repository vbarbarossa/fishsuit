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