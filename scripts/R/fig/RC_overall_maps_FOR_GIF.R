#For guidance, Nature's standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). 
#The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).

source('R/4targets/MASTER.R');

library(raster); library(foreach); library(sf); library(dplyr); library(matrixStats); library(ggplot2)
#> FUNCTIONS AND BASE LAYERS -------------------------------------------------------------------------------------

# attribute table with year/warming target/scenario/climate model combos
attr.tab <- foreach(i = seq_along(warming_targets),.combine = 'rbind') %do% {
  data.frame(
    clmod = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[1]),
    scen = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[2]),
    year = warming_tab[,i+1],
    warmt = warming_targets[i])
}
attr.tab <- attr.tab[!is.na(attr.tab$year),]
row.names(attr.tab) <- NULL

# function that takes in input the warming target and variable and gives the rasterized relative losses per grid cell
rasterize_rel_losses <- function(wtg,var = 'all',dir_root = 'proc/',dispersal=FALSE){
  
  template <- temp %>% dplyr::select(row_no,geometry)
  atab <- attr.tab[attr.tab$warmt == wtg,]
  
  tab <- foreach(x = 1:nrow(atab),.combine='cbind') %do%{
    
    if(!dispersal) t <- readRDS(paste0(dir_root,'/',atab$clmod[x],'/SR_tab_',
                                       atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
    if(dispersal) t <- readRDS(paste0(dir_root,'/',atab$clmod[x],'/SR_tab_dispersal2_',
                                      atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
    
    if(!disperal) t$val <- (t$occ - t[,var])/t$occ
    if(dispersal) t$val <- (t[,paste0('occ_',var)] - t[,var])/t[,paste0('occ_',var)]
    
    return(t[,'val'])
  }
  template$val <- rowMedians(as.matrix(tab))
  template <- template[!is.na(template$val),]
  
  return(
    raster(as(as_Spatial(template[,'val']), "SpatialPixelsDataFrame"))
  )
  
}

# base layers
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)

#> OVERALL MAPS -------------------------------------------------------------------------------------------------------
# compile raster layers
lst_ras <- list()
lst_ras_dsp <- list()
for(i in seq_along(warming_targets)){
  
  lst_ras[[i]] <- raster(paste0(dir_mod,'figshare/local_cumulative_range_loss_no_dispersal_',warming_targets[i],'.tif'))
  lst_ras_dsp[[i]] <- raster(paste0(dir_mod,'figshare/local_cumulative_range_loss_max_dispersal_',warming_targets[i],'.tif'))
  
}


# convert to ggplot format
deg <- c('1.5°C','2.0°C','3.2°C','4.5°C')
df <- rbind(
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(lst_ras[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i], scenario = 'No dispersal') %>%
      dplyr::select(x,y,value = colnames(.)[1],warming,scenario)
  },
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(lst_ras_dsp[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i], scenario = 'Maximal dispersal') %>%
      dplyr::select(x,y,value = colnames(.)[1],warming,scenario)
  }
)

df$warming <- factor(df$warming, levels = deg)
df$scenario <- factor(df$scenario, levels = c('No dispersal','Maximal dispersal'))

# and draw
for(i in seq_along(warming_targets)){
  
  p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_tile(data = df %>% filter(warming == deg[i], scenario == 'No dispersal'), aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_viridis_c(breaks = seq(0,1,0.1),
                       labels = seq(0,1,0.1),
                       limits = c(0,1),
                       option = 'C',direction = -1,na.value = "transparent") +
  # facet_grid(.~warming) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_line(color=NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 18),
        legend.title = element_blank()
  )

# ggsave(paste0(dir_mod,'figs/presentation/maps_RC_',i,'.jpg'),p,
#        width = 180,height = 100,dpi = 600,units = 'mm')
ggsave(paste0(dir_mod,'figs/presentation/maps_RC_big_',i,'.jpg'),p,
       width = 400,height = 220,dpi = 600,units = 'mm')

}
