source('R/4targets/MASTER.R')

# Library
library(fmsb); library(foreach); library(grid); library(RColorBrewer); library(dplyr); library(tidyr); library(rfishbase); library(dplyr)

#------------------------------------------------------------------------------------------------------------------------------------------------------
#> DATA

ESH_filename <- 'tabs/ESH_merged_warmingtargets_confidencemargin.rds'
if(file.exists(ESH_filename)){
  tab_long <- readRDS(ESH_filename)
}else{
  
  tab_long <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    t0 <- read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab.csv'))
    t <- cbind(t0,as.data.frame(do.call('rbind',strsplit(as.character(t0[,'comboscen']),'_'))))
    t$GCM <- clmod
    t$comboscen <- paste0(t$comboscen,'_',t$GCM)
    colnames(t)[colnames(t) %in% paste0('V',1:3)] <- c('RCP','warmt','year')
    
    return(t)
  }
  tab_long <- droplevels(tab_long)
  tab_long$GCM <- factor(tab_long$GCM)
  tab_long$comboscen <- factor(tab_long$comboscen)
  saveRDS(tab_long,ESH_filename)
}

# retrieve IUCN metadata
iucn_simplified <- foreach(i = 1:2,.combine = 'rbind') %do% foreign::read.dbf(paste0('data/FW_FISH_PART_',i,'.dbf')) %>%
  as_tibble() %>%
  distinct(binomial,.keep_all = T)

# habitat metadata
hab_data <- read.csv('data/iucn_habitat_type.csv') %>%
  as_tibble() %>%
  distinct() %>%
  select(-c(paste0('type_',1:71)))

# fishbase metadata
# set fishbase version
fishbase <- read.csv('data/iucn_fishbase.csv') %>%
  as_tibble() %>%
  mutate(binomial = iucn_name)

#hydro-climatic spatial features metadata
hydro <- read.csv('data/species_climate_feow_waterhed_area.csv')


tab <- tab_long %>%
  group_by(id_no) %>%
  summarize(
    esh1 = mean(ESH_all[warmt == '1.5C']),
    esh2 = mean(ESH_all[warmt == '2.0C']),
    esh3 = mean(ESH_all[warmt == '3.2C']),
    esh4 = mean(ESH_all[warmt == '4.5C'])
  ) %>%
  as_tibble() %>%
  inner_join(iucn_simplified,by='id_no') %>%
  left_join(hab_data,by='binomial') %>%
  left_join(hydro,by='id_no') %>%
  left_join(fishbase,by='binomial') %>%
  droplevels(.) %>%
  distinct(id_no,.keep_all = T)

#---------------------------------------------------------------------------------------------------------------------------------------
#> CATEGORIES ADJUSTMENT

# IUCN code---
tab$code[tab$code == 'LR/lc'] <- 'LC'
tab$code[tab$code == 'EX'] <- 'DD'
# levels(tab$code) <- c('Critically Endangered','Data Deficient','Endangered',
#                       'Extinct','Least Concern','Near Threatened','Vulnerable','LR/lc')

# IUCN habitat type---
tab$habtype <- 'Le'
tab$habtype[tab$lotic == 1 & tab$lentic == 0] <- 'Lo'
tab$habtype[tab$lotic == 1 & tab$lentic == 1] <- 'LL'
tab$habtype[tab$marine == 't'] <- 'M' # consider any fish that is also marine as diadromous

# koppen-geiger climatic zoning---
# plot(raster('data/Koeppen-Geiger-Classification-Reclassfied_3min_major.nc'))
# correct for 24 species with same of cells on 2 zones
tab$kg_main[tab$kg_main == '1;2'] <- 1
tab$kg_main[tab$kg_main == '1;3'] <- 1
tab <- droplevels(tab)
# assign names to levels
levels(tab$kg_main) <- LETTERS[1:5]
tab$kg_main <- factor(tab$kg_main,levels = rev(LETTERS[1:5]))
#A: Equatorial
#B: Arid
#C: Warm-temperate
#D: Snow
#E: Polar

#Migratory---
#Anacat #<<<<<<<<< there are about 52 oceanodromous fishes (spend entire life in ocean, weird!!)
levels(tab$AnaCat) <- c(NA,rep('Diad.',5),'Non.',NA,'Pota.')

#Importance---
levels(tab$Importance) <- c(NA,'Com.','HCom.','MCom.','NoInt.','NoInt.','Subs.')

#PriceCateg---
levels(tab$PriceCateg) <- c('hi','lo','med',NA,'vhi')

#FoodTroph
tab$FoodTroph[is.na(tab$FoodTroph)] <- 0
tab$foodtrophcat[tab$FoodTroph > 0] <- 'Herbi.'
tab$foodtrophcat[tab$FoodTroph > 2.19 & tab$FoodTroph <= 2.79] <- 'Omni.'
tab$foodtrophcat[tab$FoodTroph > 2.79] <- 'Carni.'
tab$FoodTroph[tab$FoodTroph == 0] <- NA

#Vulnerability #interesting to compare with the IUCN categories, e.g. boxplots
#Anacat #train an ANN to validate migratory categories inferred from IUCN with fishbase data, maybe phylogenetic approach?

tab <- droplevels(tab)

dir_out <- dir_('stats_on_traits_output/')

df <- tab  %>%
  # select covariates
  select(binomial, order = order_,family, climate = kg_main,area,length = Length,habitat = habtype,code,foodtrophcat,
         FoodTroph,importance = Importance, RC1.5 = esh1, RC2.0 = esh2, RC3.2 = esh3, RC4.5 = esh4)

write.csv(df,paste0(dir_out,'df_for_stats_on_traits.csv'),row.names = F)





#> DISPERSAL TAB--------------------------------------------------------------------------------------------------------------

ESH_filename <- 'tabs/ESH_merged_warmingtargets_confidencemargin_dispersal.rds'
if(file.exists(ESH_filename)){
  tab_long <- readRDS(ESH_filename)
}else{
  
  tab_long <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    t0 <- read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab_dispersal3.csv'))
    t <- cbind(t0,as.data.frame(do.call('rbind',strsplit(as.character(t0[,'comboscen']),'_'))))
    t$GCM <- clmod
    t$comboscen <- paste0(t$comboscen,'_',t$GCM)
    colnames(t)[colnames(t) %in% paste0('V',1:3)] <- c('RCP','warmt','year')
    
    return(t)
  }
  tab_long <- droplevels(tab_long)
  tab_long$GCM <- factor(tab_long$GCM)
  tab_long$comboscen <- factor(tab_long$comboscen)
  saveRDS(tab_long,ESH_filename)
}

# retrieve IUCN metadata
iucn_simplified <- foreach(i = 1:2,.combine = 'rbind') %do% foreign::read.dbf(paste0('data/FW_FISH_PART_',i,'.dbf')) %>%
  as_tibble() %>%
  distinct(binomial,.keep_all = T)

# habitat metadata
hab_data <- read.csv('data/iucn_habitat_type.csv') %>%
  as_tibble() %>%
  distinct() %>%
  select(-c(paste0('type_',1:71)))

# fishbase metadata
# set fishbase version
fishbase <- read.csv('data/iucn_fishbase.csv') %>%
  as_tibble() %>%
  mutate(binomial = iucn_name)

#hydro-climatic spatial features metadata
hydro <- read.csv('data/species_climate_feow_waterhed_area.csv')


tab <- tab_long %>%
  group_by(id_no) %>%
  summarize(
    esh1 = mean(ESH_all[warmt == '1.5C']),
    esh2 = mean(ESH_all[warmt == '2.0C']),
    esh3 = mean(ESH_all[warmt == '3.2C']),
    esh4 = mean(ESH_all[warmt == '4.5C'])
  ) %>%
  as_tibble() %>%
  inner_join(iucn_simplified,by='id_no') %>%
  left_join(hab_data,by='binomial') %>%
  left_join(hydro,by='id_no') %>%
  left_join(fishbase,by='binomial') %>%
  droplevels(.) %>%
  distinct(id_no,.keep_all = T)

#---------------------------------------------------------------------------------------------------------------------------------------
#> CATEGORIES ADJUSTMENT

# IUCN code---
tab$code[tab$code == 'LR/lc'] <- 'LC'
tab$code[tab$code == 'EX'] <- 'DD'
# levels(tab$code) <- c('Critically Endangered','Data Deficient','Endangered',
#                       'Extinct','Least Concern','Near Threatened','Vulnerable','LR/lc')

# IUCN habitat type---
tab$habtype <- 'Le'
tab$habtype[tab$lotic == 1 & tab$lentic == 0] <- 'Lo'
tab$habtype[tab$lotic == 1 & tab$lentic == 1] <- 'LL'
tab$habtype[tab$marine == 't'] <- 'M' # consider any fish that is also marine as diadromous

# koppen-geiger climatic zoning---
# plot(raster('data/Koeppen-Geiger-Classification-Reclassfied_3min_major.nc'))
# correct for 24 species with same of cells on 2 zones
tab$kg_main[tab$kg_main == '1;2'] <- 1
tab$kg_main[tab$kg_main == '1;3'] <- 1
tab <- droplevels(tab)
# assign names to levels
levels(tab$kg_main) <- LETTERS[1:5]
tab$kg_main <- factor(tab$kg_main,levels = rev(LETTERS[1:5]))
#A: Equatorial
#B: Arid
#C: Warm-temperate
#D: Snow
#E: Polar

#Migratory---
#Anacat #<<<<<<<<< there are about 52 oceanodromous fishes (spend entire life in ocean, weird!!)
levels(tab$AnaCat) <- c(NA,rep('Diad.',5),'Non.',NA,'Pota.')

#Importance---
levels(tab$Importance) <- c(NA,'Com.','HCom.','MCom.','NoInt.','NoInt.','Subs.')

#PriceCateg---
levels(tab$PriceCateg) <- c('hi','lo','med',NA,'vhi')

#FoodTroph
tab$FoodTroph[is.na(tab$FoodTroph)] <- 0
tab$foodtrophcat[tab$FoodTroph > 0] <- 'Herbi.'
tab$foodtrophcat[tab$FoodTroph > 2.19 & tab$FoodTroph <= 2.79] <- 'Omni.'
tab$foodtrophcat[tab$FoodTroph > 2.79] <- 'Carni.'
tab$FoodTroph[tab$FoodTroph == 0] <- NA

#Vulnerability #interesting to compare with the IUCN categories, e.g. boxplots
#Anacat #train an ANN to validate migratory categories inferred from IUCN with fishbase data, maybe phylogenetic approach?

tab <- droplevels(tab)

dir_out <- dir_('stats_on_traits_output/')

df <- tab  %>%
  # select covariates
  select(binomial, order = order_,family, climate = kg_main,area,length = Length,habitat = habtype,code,foodtrophcat,
         FoodTroph,importance = Importance, RC1.5 = esh1, RC2.0 = esh2, RC3.2 = esh3, RC4.5 = esh4)

write.csv(df,paste0(dir_out,'df_for_stats_on_traits_dispersal.csv'),row.names = F)

