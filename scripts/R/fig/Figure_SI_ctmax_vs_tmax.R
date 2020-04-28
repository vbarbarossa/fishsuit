
source('R/4targets/MASTER.R')
library(foreach); library(tidyr)

# load iucn species
iucn <- foreach(i = 1:2,.combine='rbind') %do% foreign::read.dbf(paste0('data/FW_FISH_PART_',i,'.dbf'))
iucn_simple <- iucn[,c('id_no','binomial')][!duplicated(iucn[,c('id_no','binomial')]),]
row.names(iucn_simple) <- 1:nrow(iucn_simple)

#> RETRIEVE CTMAX LAB DATA --------------------------------------------------------------------------------------

# based on observed data
tab <- read.csv('data/Comte__Olden_original.csv') %>% .[.$Realm.affinity == 'Freshwater',]
length(unique(tab$Species)) # 327 as in the paper

# try with a simple average over the values available for each species
lab_data <- foreach(sp = unique(as.character(tab$Species)),.combine = 'rbind') %do% {
  tt <- tab[tab$Species == sp,]
  data.frame(
    Species = sp,
    Acclimation_lab = mean(tt$Temperature.of.acclimation...C.,na.rm=T),
    Acclimation.SD_lab = sd(tt$Temperature.of.acclimation...C.,na.rm=T),
    CTmax_lab = mean(tt$Thermal.limit...C.,na.rm=T),
    CTmax_max = max(tt$Thermal.limit...C.,na.rm=T),
    CTmax_min = min(tt$Thermal.limit...C.,na.rm=T),
    CTmax_sd =  sd(tt$Thermal.limit...C.,na.rm=T),
    CTmax.SD_lab = sd(tt$Thermal.limit...C.,na.rm=T)
    
  )
}

#> RETRIEVE TMAX MODELLED DATA --------------------------------------------------------------------------------------

climate_models <- c('gfdl','ipsl','hadgem','miroc','noresm') #,'noresm'

# plot against Tmax data
dat <- merge(iucn_simple,lab_data,by.x='binomial',by.y = 'Species') # 214

data <- dat[dat$id_no %in% read.csv(paste0(dir_mod,'/miroc/niches/niches_filtered.csv'))[,'IUCN_ID'] ,]
# 1 lost

# retrieve 97.5% of Tmax from different climate models
tc <- cbind(sort(readRDS(paste0(dir_mod,'/miroc/niches/Tma.rds'))[,c('IUCN_ID')]),
            foreach(clmod = climate_models,.combine = 'cbind') %do% {
              
              t <- readRDS(paste0(dir_mod,clmod,'/niches/Tma.rds'))[,c('IUCN_ID','97.5%')]
              t[,2] <- t[,2:ncol(t)] - 273.15
              t <- t[order(t$IUCN_ID),]
              
              return(data.frame(t[,2]))
              
            })
colnames(tc) <- c('IUCN_ID',climate_models)
tc$Tmax_median <- apply(tc[2:(length(climate_models)+1)],1,median)
tc$Tmax_mean <- apply(tc[2:(length(climate_models)+1)],1,mean)
tc$Tmax_max <- apply(tc[2:(length(climate_models)+1)],1,max)
tc$Tmax_min <- apply(tc[2:(length(climate_models)+1)],1,min)

data <- merge(data,tc,by.x = 'id_no',by.y='IUCN_ID')
mse_val <- MSE(data$CTmax_lab,data$Tmax_median)
rsq_val <- r.squared(data$CTmax_lab,data$Tmax_median)

esh_merged_warming_targets <- readRDS('tabs/ESH_merged_warmingtargets.rds')
addt <- foreach(i = data$id_no,.combine = 'rbind') %do% {
  t <- esh_merged_warming_targets[,c('id_no','no.cells')]
  t[t$id_no == i,][1,]
}
data <- merge(data,addt,by='id_no')

# merge habitat data
hab_type <- read.csv('data/iucn_habitat_type.csv')
hab_data <- foreach(i = as.character(unique(data$binomial)),.combine = 'rbind') %do% hab_type[hab_type$binomial == i,][1,] #lookslike there are more rows than species

data <- merge(data,hab_data[,c('binomial','lotic','lentic','seasonal')],by='binomial') # we lose 76 species when merging <<<<<<

data$hab_type <- 'Lentic'
data$hab_type[data$lotic == 1] <- 'Lotic' # consider any fish that lives in rivers and not lakes as potadromous
data$hab_type <- factor(data$hab_type,levels = c('Lotic','Lentic'))


library(ggplot2)
p <- ggplot(data,aes(x=CTmax_lab,y=Tmax_median)) +
  geom_abline(color = 'black',linetype='dashed') +
  geom_smooth(method='lm',formula=y~x,fullrange = T,linetype = 'dashed') +
  geom_errorbar(aes(ymin = Tmax_min, ymax = Tmax_max), color = 'gray',linetype='solid',width = 0.3,alpha = 0.7) +
  geom_errorbarh(aes(xmin = CTmax_min, xmax = CTmax_max), color = 'gray',linetype='solid',height = 0.3,alpha = 0.7) +
  geom_point(aes(size = log(no.cells,base = 10),shape = hab_type),alpha=0.5) +
  xlab(expression('Critical Thermal Maximum from lab data '*'['^o*C*']')) +
  ylab(expression('Upper threshold of maximum weekly water temperature '*'['^o*C*']')) +
  coord_cartesian(xlim = c(20,45),ylim = c(20,45),expand=T) +
  guides(size = FALSE,shape = FALSE) +
  theme_bw() +
  theme(
    text = element_text(size = 15,color='black'),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    panel.border = element_rect(color='black')
  )
p

lm1 <- lm(Tmax_median ~ CTmax_lab, data = data)
lm2 <- lm(Tmax_median ~ 0 + CTmax_lab, data = data)
anova(lm1,lm2)

ggsave(paste0(dir_mod,'figs/Figure_SI_ctmax_vs_tmax_herrbar.jpeg'),p,
       width = 120,height = 120,units = 'mm',dpi = 300,scale = 1.5 )

# no trendline and climate
tab <- read.csv('data/species_climate_feow_waterhed_area.csv')
tab$kg_main[tab$kg_main == '1;2'] <- 1
tab$kg_main[tab$kg_main == '1;3'] <- 1
tab <- droplevels(tab)

df <- droplevels(merge(data,tab,by='id_no'))

levels(df$kg_main) <- LETTERS[1:5]


library(ggplot2)
p <- ggplot(df,aes(x=CTmax_lab,y=Tmax_median)) +
  geom_abline(color = 'black',linetype='dashed') +
  geom_errorbar(aes(ymin = Tmax_min, ymax = Tmax_max), color = 'gray',linetype='solid',width = 0.3,alpha = 0.7) +
  geom_point(aes(size = log(area,base = 10),shape = hab_type,color = kg_main),alpha=0.5) +
  xlab(expression('Critical Thermal Maximum from lab data '*'['^o*C*']')) +
  ylab(expression('Upper threshold of maximum weekly water temperature '*'['^o*C*']')) +
  coord_cartesian(xlim = c(20,45),ylim = c(20,45),expand=T) +
  guides(size = guide_legend(title="Range\n[log10-km2]"),shape = guide_legend(title="Habitat"),color = guide_legend(title="Climate")) +
  theme_bw() +
  theme(
    text = element_text(size = 15,color='black'),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    panel.border = element_rect(color='black')
  )
p

ggsave(paste0(dir_mod,'figs/Figure_SI_ctmax_vs_tmax_categories.jpeg'),p,
       width = 140,height = 120,units = 'mm',dpi = 600,scale = 1.5 )

#----------------------------------------------------------------------
# simple for presentation
library(ggplot2)
p <- ggplot(df,aes(x=CTmax_lab,y=Tmax_median)) +
  geom_abline(color = 'black',linetype='dashed') +
  geom_errorbar(aes(ymin = Tmax_min, ymax = Tmax_max), color = 'gray',linetype='solid',width = 0.3,alpha = 0.7) +
  geom_point(alpha=0.5,size=5) +
  xlab(expression('Lab data '*'['^o*C*']')) +
  ylab(expression('This study '*'['^o*C*']')) +
  coord_cartesian(xlim = c(20,45),ylim = c(20,45),expand=T) +
  theme_bw() +
  theme(
    text = element_text(size = 20,color='black'),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    panel.border = element_rect(color='black')
  )
p

ggsave(paste0(dir_mod,'figs/presentation/ctmax_vs_tmax.jpeg'),p,
       width = 120,height = 120,units = 'mm',dpi = 600,scale = 1.5 )
