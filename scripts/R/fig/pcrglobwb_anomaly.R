source('config.R');

source_layers <- 'proc/pcrglobwb_anomaly_ensemble_perc_change.rds'

if(file.exists(source_layers)){
  diff_ens <- readRDS(source_layers)
}else{
  
  library(raster); library(dplyr)
  
  # attribute table with year/warming target/scenario/climate model combos
  library(foreach)
  attr.tab <- foreach(i = seq_along(warming_targets),.combine = 'rbind') %do% {
    data.frame(
      clmod = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[1]),
      scen = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[2]),
      year = warming_tab[,i+1],
      warmt = warming_targets[i])
  }
  attr.tab <- attr.tab[!is.na(attr.tab$year),]
  row.names(attr.tab) <- NULL
  
  # enseble historical output
  ensemble_median <- function(ls.files){
    library(foreach); library(dplyr); library(raster)
    r <- foreach(i = ls.files) %do% raster(i) %>% brick()
    
    beginCluster(8)
    med <- clusterR(r, calc, args=list(median, na.rm=T))
    endCluster()
    
    return(med)
    
  }
  
  vars <- c('Qmi','Qzf','Tma')
  
  hist <- lapply(
    foreach(var = vars) %do% paste0('proc/',climate_models,'/pcrglobwb_processed/merged/',var,'_hist.tif'),
    ensemble_median
  )
  names(hist) <- vars
  
  diff_ens <- foreach(wt = unique(attr.tab$warmt)) %do% {
    
    attrs <- droplevels(attr.tab[attr.tab$warmt == wt,])
    
    fut <- lapply(
      foreach(var = vars) %do% {
        foreach(r = 1:nrow(attrs),.combine = 'c') %do% paste0('proc/',attrs[r,'clmod'],
                                                              '/pcrglobwb_processed/merged/',var,'_',
                                                              attrs[r,'scen'],'_',attrs[r,'warmt'],'C_',
                                                              attrs[r,'year'],'.tif')
      },
      ensemble_median
    )
    names(fut) <- vars
    
    diff <- lapply(1:3,function(x) {(mask(fut[[x]],hist[[x]] >= 0) - mask(hist[[x]],hist[[x]] >= 0))/mask(hist[[x]],hist[[x]] >= 0)})
    names(diff) <- vars
    
    return(diff)
    
  }
  
  names(diff_ens) <- warming_targets
  saveRDS(diff_ens,source_layers)
  
}

#--------------------------------------------------------------------------------------
# one plot for each variable at 4 different warming targets
library(ggplot2); library(RColorBrewer); library(sf)

# base layers
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)


# clean raster layers?


# convert to ggplot format
deg <- c('1.5°C','2.0°C','3.2°C','4.5°C')
df <- bind_rows(
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(diff_ens[[warming_targets[i]]][['Qzf']] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i], var = 'Qzf',varlong = 'max. no. zero-flow weeks') %>%
      dplyr::select(x,y,value = layer,warming,var,varlong) %>%
      filter(!is.na(value))
  } %>% as_tibble(),
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(diff_ens[[warming_targets[i]]][['Qmi']] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i], var = 'Qmi',varlong = 'min. weekly flow') %>%
      dplyr::select(x,y,value = layer,warming,var,varlong) %>%
      filter(!is.na(value))
  } %>% as_tibble(),
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(diff_ens[[warming_targets[i]]][['Tma']] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i], var = 'Tma',varlong = 'max. weekly water temp.') %>%
      dplyr::select(x,y,value = layer,warming,var,varlong) %>%
      filter(!is.na(value))
  } %>% as_tibble()
)

df$warming <- factor(df$warming, levels = deg)
df$var <- factor(df$var, levels = c('Qmi','Qzf','Tma'))
df$varlong <- factor(df$varlong, levels = c('min. weekly flow',
                                            'max. no. zero-flow weeks',
                                            'max. weekly water temp.'))

# do one plot per variable to manage better the colorscale
# water temperature has much smaller variations in percentage of change than flow

# QZF--------------------------------------------------------------------------
dd <- df %>% filter(var == 'Qzf') %>% mutate(value = round(value*100,0))
dd$value[dd$value == Inf] <- max(dd$value[!is.infinite(dd$value)])
dd$value[dd$value == -Inf] <- min(dd$value[!is.infinite(dd$value)])

dd$value[dd$value > 100] <- 100

# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = 'grey90', lwd = NA) + # '#e4cead' yellow ochre
  geom_tile(data = dd, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_gradient2(mid = 'grey90',
                       midpoint = 0,
                       low = scales::muted('blue'),
                       high = scales::muted('red'),
                       breaks = seq(-100,100,25),
                       labels = paste0(c(seq(-100,0,25),paste0('+',seq(25,100,25))),'%'),
                       na.value = 'transparent') +
  facet_wrap('warming',nrow = 2) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        panel.spacing = unit(0.1, "lines"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.title = element_blank()
  )
# p

ggsave(paste0(dir_mod,'figs/maps_pcrglobwb_anomaly_Qzf.jpg'),p,
       width = 200,height = 130,dpi = 600,units = 'mm')

# QMIN--------------------------------------------------------------------------
dd <- df %>% filter(var == 'Qmi') %>% mutate(value = round(value*100,0))
dd$value[dd$value == Inf] <- max(dd$value[!is.infinite(dd$value)])
dd$value[dd$value == -Inf] <- min(dd$value[!is.infinite(dd$value)])

dd$value[dd$value > 100] <- 100


# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = 'grey90', lwd = NA) + # '#e4cead' yellow ochre
  geom_tile(data = dd, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_gradient2(mid = 'grey90',
                       midpoint = 0,
                       breaks = seq(-100,100,25),
                       labels = paste0(c(seq(-100,0,25),paste0('+',seq(25,100,25))),'%'),
                       na.value = 'transparent') +
  facet_wrap('warming',nrow = 2) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        panel.spacing = unit(0.1, "lines"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.title = element_blank()
  )
# p

ggsave(paste0(dir_mod,'figs/maps_pcrglobwb_anomaly_Qmi.jpg'),p,
       width = 200,height = 130,dpi = 600,units = 'mm')



# TMAX-----------------------------------------------------------------------
dd <- df %>% filter(var == 'Tma') %>% mutate(value = round(value*100,0))
dd$value[dd$value == Inf] <- max(dd$value[!is.infinite(dd$value)])
dd$value[dd$value == -Inf] <- min(dd$value[!is.infinite(dd$value)])

dd$value[dd$value > 5] <- 6
dd$value[dd$value < -5] <- -6


# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = 'grey90', lwd = NA) + # '#e4cead' yellow ochre
  geom_tile(data = dd, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_gradient2(mid = 'grey90',
                       low = scales::muted('blue'),
                       high = scales::muted('red'),
                       midpoint = 0,
                       breaks = seq(-6,6,1),
                       labels = paste0(c('<5',seq(-5,0,1),paste0('+',seq(1,5,1)),'>5'),'%'),
                       na.value = 'transparent') +
  facet_wrap('warming',nrow = 2) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        panel.spacing = unit(0.1, "lines"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.title = element_blank()
  )
# p

ggsave(paste0(dir_mod,'figs/maps_pcrglobwb_anomaly_Tma.jpg'),p,
       width = 200,height = 130,dpi = 600,units = 'mm')

