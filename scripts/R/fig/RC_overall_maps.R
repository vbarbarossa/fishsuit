#For guidance, Nature's standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). 
#The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).

source('config.R'); # always load as functions are loaded within this script

library(raster); library(foreach); library(sf); library(dplyr); library(matrixStats); library(ggplot2)
#> FUNCTIONS AND BASE LAYERS -------------------------------------------------------------------------------------

# template points
temp <- readRDS('proc/ssp/points_template.rds')

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

# if already compiled, use the shortcut:
# lst_ras <- list()
# lst_ras_dsp <- list()
# for(i in seq_along(warming_targets)){
#   
#   lst_ras[[i]] <- raster(paste0('figshare/local_cumulative_range_loss_no_dispersal_',warming_targets[i],'.tif'))
#   lst_ras_dsp[[i]] <- raster(paste0('figshare/local_cumulative_range_loss_max_dispersal_',warming_targets[i],'.tif'))
#   
# }
# names(lst_ras) <- warming_targets
# names(lst_ras_dsp) <- warming_targets


lst_ras <- lapply(warming_targets,function(x) rasterize_rel_losses(wtg = x))
names(lst_ras) <- warming_targets

lst_ras_dsp <- lapply(warming_targets,function(x) rasterize_rel_losses(wtg = x,dispersal = TRUE))
names(lst_ras_dsp) <- warming_targets

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
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_raster(data = df, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_viridis_c(breaks = seq(0,1,0.1),
                       labels = seq(0,1,0.1),
                       limits = c(0,1),
                       option = 'C',na.value = "transparent") +
  facet_grid(warming~scenario) +
  theme_minimal() +
  theme(text = element_text(size = 13),
        panel.grid.major = element_line(color=NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 13),
        legend.title = element_blank()
  )
# p
ggsave(paste0('figs/maps_RC.jpg'),p,
       width = 183,height = 200,dpi = 600,units = 'mm')
# ggsave(paste0('figs/maps_RC.pdf'),p,
#        width = 183,height = 200,units = 'mm')

# save source files for figshare
for(i in seq_along(warming_targets)){
  
  writeRaster(lst_ras[[i]],paste0('figshare/local_cumulative_range_loss_no_dispersal_',warming_targets[i],'.tif'),format = 'GTiff',overwrite = T)
  writeRaster(lst_ras_dsp[[i]],paste0('figshare/local_cumulative_range_loss_max_dispersal_',warming_targets[i],'.tif'),format = 'GTiff',overwrite = T)
  
}

# no dispersal only
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_raster(data = df %>% filter(scenario == 'No dispersal'), 
              aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_viridis_c(breaks = seq(0,1,0.1),
                       labels = seq(0,1,0.1),
                       limits = c(0,1),
                       option = 'C',na.value = "transparent") +
  facet_wrap('warming',nrow = 2) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(5,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.title = element_blank()
  )
# p
ggsave(paste0('figs/maps_RC_nodisp.jpg'),p,
       width = 183,height = 120,dpi = 600,units = 'mm')
# ggsave(paste0('figs/maps_RC.pdf'),p,
#        width = 183,height = 200,units = 'mm')


#> CONTRIBUTION T & Q --------------------------------------------------------------------------------------

for(w in seq_along(warming_targets)){
  # compile raster layers, list of variables for 3.2deg
  varnames <- c('Q_all','T_all','both_QT')
  lst_ras <- lapply(varnames,function(x) rasterize_rel_losses(wtg = warming_targets[w],var = x))
  names(lst_ras) <- varnames
  
  lst_ras_dsp <- lapply(varnames,function(x) rasterize_rel_losses(wtg = warming_targets[w],var = x,dispersal = TRUE))
  names(lst_ras_dsp) <- varnames
  
  # convert to ggplot format
  df <- rbind(
    foreach(i = seq_along(varnames),.combine='rbind') %do% {
      as(lst_ras[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
        as.data.frame(.) %>%
        mutate(varname = varnames[i], scenario = 'No dispersal') %>%
        dplyr::select(x,y,value = val,varname,scenario)
    },
    foreach(i = seq_along(varnames),.combine='rbind') %do% {
      as(lst_ras_dsp[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
        as.data.frame(.) %>%
        mutate(varname = varnames[i], scenario = 'Maximal dispersal') %>%
        dplyr::select(x,y,value = val,varname,scenario)
    }
  ) %>%
    mutate(varname = forcats::fct_recode(varname, 'Streamflow' = 'Q_all', 'Water temperature' = 'T_all','Both' = 'both_QT')) %>%
    mutate(varname = factor(varname, levels = c('Water temperature', 'Streamflow', 'Both')),
           scenario = factor(scenario, levels = c('No dispersal','Maximal dispersal')))
  
  # and draw
  p <- ggplot() +
    geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
    geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
    geom_sf(data = world, fill = "grey90", lwd = NA) +
    geom_raster(data = df, aes(x=x, y=y, fill=value), alpha=0.8) +
    scale_fill_viridis_c(breaks = seq(0,1,0.1),
                         labels = seq(0,1,0.1),
                         limits = c(0,1),
                         option = 'C',na.value = "transparent") +
    facet_grid(varname~scenario) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          panel.grid.major = element_line(color=NA),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = 'bottom',
          legend.key.width = unit(6,'line'),
          strip.background = element_rect('white'),
          strip.background.x = element_blank(),
          strip.background.y = element_blank(),
          strip.text = element_text(angle = 0, vjust = -1, size = 13),
          legend.title = element_blank()
    )
  # p
  ggsave(paste0('figs/maps_RC_QTcontribution_',warming_targets[w],'.jpg'),p,
         width = 183,height = 170,dpi = 600,units = 'mm')
  # ggsave(paste0('figs/maps_RC.pdf'),p,
  #        width = 183,height = 200,units = 'mm')
  
  # save source files for figshare
  for(i in seq_along(varnames)){
    
    writeRaster(lst_ras[[i]],paste0('figshare/local_cumulative_range_loss_no_dispersal_',warming_targets[w],'_',varnames[i],'.tif'),format = 'GTiff',overwrite = T)
    writeRaster(lst_ras_dsp[[i]],paste0('figshare/local_cumulative_range_loss_max_dispersal_',warming_targets[w],'_',varnames[i],'.tif'),format = 'GTiff',overwrite = T)
    
  }
  
  
}


#> CONTRIBUTION single variables --------------------------------------------------------------------------------------

for(w in seq_along(warming_targets)){
  
  # compile raster layers, list of variables for 3.2deg
  varnames <- c("Qmi","Qzf","Qma","Tma","Tmi")
  
  lst_ras <- lapply(varnames,function(x) rasterize_rel_losses(wtg = warming_targets[w],var = x))
  names(lst_ras) <- varnames
  
  lst_ras_dsp <- lapply(varnames,function(x) rasterize_rel_losses(wtg = warming_targets[w],var = x,dispersal = TRUE))
  names(lst_ras_dsp) <- varnames
  
  # convert to ggplot format
  df <- rbind(
    foreach(i = seq_along(varnames),.combine='rbind') %do% {
      as(lst_ras[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
        as.data.frame(.) %>%
        mutate(varname = varnames[i], scenario = 'No dispersal') %>%
        dplyr::select(x,y,value = val,varname,scenario)
    },
    foreach(i = seq_along(varnames),.combine='rbind') %do% {
      as(lst_ras_dsp[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
        as.data.frame(.) %>%
        mutate(varname = varnames[i], scenario = 'Maximal dispersal') %>%
        dplyr::select(x,y,value = val,varname,scenario)
    }
  ) %>%
    # mutate(varname = forcats::fct_recode(varname, 'Streamflow' = 'Q_all', 'Water temperature' = 'T_all','Both' = 'both_QT')) %>%
    mutate(scenario = factor(scenario, levels = c('No dispersal','Maximal dispersal')))
  
  # and draw
  p <- ggplot() +
    geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
    geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
    geom_sf(data = world, fill = "grey90", lwd = NA) +
    geom_raster(data = df, aes(x=x, y=y, fill=value), alpha=0.8) +
    scale_fill_viridis_c(breaks = seq(0,1,0.1),
                         labels = seq(0,1,0.1),
                         limits = c(0,1),
                         option = 'C',na.value = "transparent") +
    facet_grid(varname~scenario) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          panel.grid.major = element_line(color=NA),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = 'bottom',
          legend.key.width = unit(6,'line'),
          strip.background = element_rect('white'),
          strip.background.x = element_blank(),
          strip.background.y = element_blank(),
          strip.text = element_text(angle = 0, vjust = -1, size = 13),
          legend.title = element_blank()
    )
  # p
  ggsave(paste0('figs/maps_RC_single_variables_',warming_targets[w],'.jpg'),p,
         width = 183,height = 170*4/3,dpi = 600,units = 'mm')
  # ggsave(paste0('figs/maps_RC.pdf'),p,
  #        width = 183,height = 200,units = 'mm')
  
  # save source files for figshare
  for(i in seq_along(varnames)){
    
    writeRaster(lst_ras[[i]],paste0('figshare/local_cumulative_range_loss_no_dispersal_',warming_targets[w],'_',varnames[i],'.tif'),format = 'GTiff',overwrite = T)
    writeRaster(lst_ras_dsp[[i]],paste0('figshare/local_cumulative_range_loss_max_dispersal_',warming_targets[w],'_',varnames[i],'.tif'),format = 'GTiff',overwrite = T)
    
  }
  
  
}


#> INITIAL RICHNESS -----------------------------------------------------------------------------------------------------

lst_ras <- foreach(source = c('hadgem/SR_tab_rcp8p5_3.2C_2058.rds',
                              'hadgem/SR_tab_dispersal_rcp8p5_3.2C_2058.rds')
) %do%{
  
  ttemp <- left_join(temp,readRDS(paste0('proc/',source)),by = 'row_no') %>% 
    dplyr::select(occ) %>% 
    filter(!is.na(occ)) 
    
  return(
    raster(as(as_Spatial(ttemp), "SpatialPixelsDataFrame"))
  )
  
}

# convert to ggplot format
deg <- c('No Dispersal','Maximal dispersal')
df <- foreach(i = seq_along(deg),.combine='rbind') %do% {
  as(lst_ras[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
    as.data.frame(.) %>%
    mutate(habitat = deg[i]) %>%
    dplyr::select(x,y,value = occ,habitat)
}

df$habitat <- factor(df$habitat, levels = deg)


for(i in deg){
  
  tt <- df
  tt$value <- tt$value %>% round(.,0) %>% log10()
  
  p <- ggplot() +
    geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
    geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
    geom_sf(data = world, fill = "grey90", lwd = NA) +
    geom_raster(data = tt[tt$habitat == i,], aes(x=x, y=y, fill=value)) +
    scale_fill_viridis_c(limits = range(tt$value),
                         na.value = "transparent",direction = -1) +
    facet_grid(habitat~.) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          panel.grid.major = element_line(color=NA),
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
  
  ggsave(paste0('figs/maps_SR_',i,'.jpg'),p,
         width = 200,height = 100,dpi = 600,units = 'mm')
  
}
for(i in seq_along(deg)){
  
  writeRaster(lst_ras[[i]],paste0('figshare/initial_species_richnees_',deg[i],'.tif'),format = 'GTiff',overwrite = T)
  
}

df$value <- df$value %>% round(.,0) %>% log10()

p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_raster(data = df, aes(x=x, y=y, fill=value)) +
  scale_fill_viridis_c(limits = range(df$value),
                       na.value = "transparent",direction = -1) +
  facet_grid(habitat~.) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
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


ggsave(paste0('figs/maps_SR.jpg'),p,
       width = 200,height = 200,dpi = 600,units = 'mm')
