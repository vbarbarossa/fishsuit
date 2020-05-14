
source('config.R')
library(foreach); library(tidyr); library(dplyr)

# load species
species <- read.csv('proc/thresholds_average_filtered.csv') %>% select(id_no,binomial,area)

#> RETRIEVE CTMAX LAB DATA --------------------------------------------------------------------------------------

# based on observed data
tab <- read.csv('data/Comte__Olden_Data_original.csv') %>% .[.$Realm.affinity == 'Freshwater',] %>% as_tibble()
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
} %>% as_tibble()

#> RETRIEVE TMAX MODELLED DATA --------------------------------------------------------------------------------------

climate_models <- c('gfdl','ipsl','hadgem','miroc','noresm') #,'noresm'

# plot against Tmax data
data <- merge(species,lab_data,by.x='binomial',by.y = 'Species') %>% as_tibble() # 293

# retrieve 97.5% of Tmax from different climate models
tc_list <- foreach(clmod = climate_models) %do% {
  
  t <- readRDS(paste0('proc/',clmod,'/niches/Tma.rds'))[,c('ID','97.5%')]
  t[,2] <- t[,2] - 273.15
  colnames(t) <- c('id_no',clmod)
  return(t)
}

# merge by id_no
tc <- tc_list[[1]] %>% as_tibble()
for(i in seq_along(tc_list)[-1]){
  tc <- inner_join(tc,tc_list[[i]])
}

# calculate stats
tc$Tmax_median <- apply(tc[2:(length(climate_models)+1)],1,median)
tc$Tmax_mean <- apply(tc[2:(length(climate_models)+1)],1,mean)
tc$Tmax_max <- apply(tc[2:(length(climate_models)+1)],1,max)
tc$Tmax_min <- apply(tc[2:(length(climate_models)+1)],1,min)

data <- inner_join(data,tc %>% mutate(id_no = as.integer(as.character(id_no))))
mse_val <- valerioUtils::MSE(data$CTmax_lab,data$Tmax_median)
rsq_val <- valerioUtils::r.squared(data$CTmax_lab,data$Tmax_median)

cor(data$CTmax_lab,data$Tmax_median)

library(ggplot2)
p <- ggplot(data,aes(x=CTmax_lab,y=Tmax_median)) +
  geom_abline(color = 'black',linetype='dashed') +
  geom_smooth(method='lm',formula=y~x,fullrange = T,linetype = 'dashed') +
  geom_errorbar(aes(ymin = Tmax_min, ymax = Tmax_max), color = 'gray',linetype='solid',width = 0.3,alpha = 0.7) +
  geom_errorbarh(aes(xmin = CTmax_min, xmax = CTmax_max), color = 'gray',linetype='solid',height = 0.3,alpha = 0.7) +
  geom_point(aes(size = log(area,base = 10)),alpha=0.5) +
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

ggsave(paste0('figs/ctmax_vs_tmax_herrbar.jpeg'),p,
       width = 120,height = 120,units = 'mm',dpi = 300,scale = 1.5 )

# no trendline and climate
tab <- read.csv('proc/species_traits.csv') %>% as_tibble()

df <- inner_join(data,tab,by = 'id_no')
df$climate_zone <- as.factor(df$climate_zone)
levels(df$climate_zone) <- c('Equatorial','Arid','Warm Temperate','Snow')
library(ggplot2)
p <- ggplot(df,aes(x=CTmax_lab,y=Tmax_median)) +
  geom_abline(color = 'black',linetype='dashed') +
  geom_errorbar(aes(ymin = Tmax_min, ymax = Tmax_max), color = 'gray',linetype='solid',width = 0.3,alpha = 0.7) +
  geom_point(aes(size = log(area.x,base = 10),color = climate_zone),alpha=0.5,shape = 16) +
  xlab(expression('Critical Thermal Maximum from lab data '*'['^o*C*']')) +
  ylab(expression('Upper threshold of maximum weekly water temperature '*'['^o*C*']')) +
  scale_color_manual(values = RColorBrewer::brewer.pal(name = 'Set1',n=9)[c(1,5,3,4)]) +
  coord_cartesian(xlim = c(20,45),ylim = c(20,45),expand=T) +
  guides(size = guide_legend(title="Range\n[log10-km2]"),color = guide_legend(title="Climate",override.aes = list(size=5))) +
  theme_bw() +
  theme(
    text = element_text(size = 15,color='black'),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    panel.border = element_rect(color='black')
  )
p

ggsave('figs/ctmax_vs_tmax.jpeg',p,
       width = 150,height = 120,units = 'mm',dpi = 600,scale = 1.5 )

#----------------------------------------------------------------------
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< not checked
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
