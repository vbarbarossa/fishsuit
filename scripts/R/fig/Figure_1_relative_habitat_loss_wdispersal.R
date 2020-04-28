#For guidance, Nature's standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). 
#The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).

source('R/4targets/MASTER.R');

library(raster); library(foreach); library(sf); library(dplyr)
#> ENSEMBLE RASTER LAYERS -------------------------------------------------------------------------------------

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

# function that takes in input the row of the attr.tab and gives the rasterized relative losses per grid cell
rasterize_rel_losses <- function(x,dir_root = 'fishsuit_completeRun_warming_4targets',dispersal=FALSE){
  if(!dispersal) tab <- readRDS(paste0(dir_root,'/',atab$clmod[x],'/SR_tab_',
                                       atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
  if(dispersal) tab <- readRDS(paste0(dir_root,'/',atab$clmod[x],'/SR_tab_',
                                      atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'_dispersal.rds'))
  
  tab <- tab[!is.na(tab@data$occ),]
  tab@data <- data.frame(rel_all = (tab@data$occ - tab@data$all)/tab@data$occ)
  return(
    rasterize(tab,raster(ext=extent(tab),resolution=1/12),field = 'rel_all')
  )
}

lst_ras_name <- paste0(dir_mod,dir_('figs/data/'),'raster_layers_figure_1.rds')

if(file.exists(lst_ras_name)){
  lst_ras <- readRDS(lst_ras_name)
}else{
  lst_ras <- foreach(n = seq_along(warming_targets)) %do% {
    
    atab <- attr.tab[attr.tab$warmt == warming_targets[n],]
    
    ncores = 7
    # parallelized ---
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    br <- brick(
      foreach(i = 1:nrow(atab),.packages='raster') %dopar% rasterize_rel_losses(x=i)
    )
    parallel::stopCluster(cl)
    # ---
    return(
      calc(br,function(x) median(x,na.rm=T))
    )
    
  }
  names(lst_ras) <- warming_targets
  
  saveRDS(file = lst_ras_name,object = lst_ras)
  
}

lst_ras_name_dsp <- paste0(dir_mod,dir_('figs/data/'),'raster_layers_dispersal_figure_1.rds')

if(file.exists(lst_ras_name_dsp)){
  lst_ras_dsp <- readRDS(lst_ras_name_dsp)
}else{
  lst_ras_dsp <- foreach(n = seq_along(warming_targets)) %do% {
    
    atab <- attr.tab[attr.tab$warmt == warming_targets[n],]
    
    ncores = 7
    # parallelized ---
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    br <- brick(
      foreach(i = 1:nrow(atab),.packages='raster') %dopar% rasterize_rel_losses(x=i,dispersal = TRUE)
    )
    parallel::stopCluster(cl)
    # ---
    return(
      calc(br,function(x) median(x,na.rm=T))
    )
    
  }
  names(lst_ras_dsp) <- warming_targets
  
  saveRDS(file = lst_ras_name_dsp,object = lst_ras_dsp)
  
}

# start ~ 2:20PM

#> MAPS -------------------------------------------------------------------------------------------------------

# base layers-----------------------------------------------------------------------------------------------------------------
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)

# convert to ggplot format--------------------------------------------------------------------------------------------------
deg <- c('1.5°C','2.0°C','3.2°C','4.5°C')
df <- rbind(
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(lst_ras[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i], scenario = 'no dispersal') %>%
      select(x,y,value = layer,warming,scenario)
  },
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(lst_ras_dsp[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i], scenario = 'unlimited dispersal') %>%
      select(x,y,value = layer,warming,scenario)
  }
)

df$warming <- factor(df$warming, levels = deg)
df$scenario <- factor(df$scenario, levels = c('no dispersal','unlimited dispersal'))

# and draw-----------------------------------------------------------------------------------------------------

library(ggplot2)
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_tile(data = df, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_viridis_c(breaks = seq(0,1,0.1),
                       labels = seq(0,1,0.1),
                       limits = c(0,1),
                       option = 'C',na.value = "transparent") +
  facet_grid(warming~scenario) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        panel.grid.major = element_line(color=NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 16),
        legend.title = element_blank()
  )
# p
ggsave(paste0(dir_mod,'figs/maps_RC.jpg'),p,
       width = 183,height = 200,dpi = 600,units = 'mm')
# ggsave(paste0(dir_mod,'figs/maps_RC.pdf'),p,
#        width = 183,height = 200,units = 'mm')


#> VIOLIN PLOTS -----------------------------------------------------------------------------------------
palette <- c(
  brewer.pal(9, "Greys")[3],
  colorRampPalette(brewer.pal(9, "YlOrRd"))(40)
)

# this is available from 'Figure_SI_boxplot_scenario.R'
tab <- readRDS(paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds'))
tab <- droplevels( tab[tab$ESH_type == 'Total',] )


tab_median <- foreach(sp = unique(tab$id_no),.combine = 'rbind') %do% {
  tsp <- droplevels(tab[tab$id_no == sp,])
  foreach(wt = paste0(warming_targets,'C'),.combine = 'rbind') %do% {
    
    t <- droplevels(tsp[tsp$warmt == wt,])
    return(
      data.frame(
        id_no = sp,
        no.cells = unique(t$no.cells),
        ESH_median = median(t$ESH,na.rm=T), #<<<<<<<<<<<<<<<<<<<< NEED TO CHECK WHY MEDIAN HERE
        warmt = wt
      )
    )
  } 
}

# this is available from 'Figure_SI_boxplot_scenario.R'
tab <- readRDS(paste0(dir_mod,'tabs/ESH_merged_warmingtargets_dispersal.rds'))
tab <- droplevels( tab[tab$ESH_type == 'Total',] )


tab_median_dsp <- foreach(sp = unique(tab$id_no),.combine = 'rbind') %do% {
  tsp <- droplevels(tab[tab$id_no == sp,])
  foreach(wt = paste0(warming_targets,'C'),.combine = 'rbind') %do% {
    
    t <- droplevels(tsp[tsp$warmt == wt,])
    return(
      data.frame(
        id_no = sp,
        no.cells = unique(t$no.cells),
        ESH_median = median(t$ESH,na.rm=T), #<<<<<<<<<<<<<<<<<<<< NEED TO CHECK WHY MEDIAN HERE
        warmt = wt
      )
    )
  } 
}

df_bp <- rbind(
  tab_median %>% mutate(scenario = 'no dispersal'),
  tab_median_dsp %>% mutate(scenario = 'dispersal')
)

df_bp$scenario <- factor(df_bp$scenario, levels = c('no dispersal','dispersal'))

library(ggplot2)
p <- ggplot(df_bp,aes(x=warmt,y=ESH_median)) + #
  geom_violin(aes(fill = scenario),lwd = .5,color='transparent') + #
  geom_boxplot(aes(color = scenario),fill='white',outlier.color = 'transparent',
               width = 0.08,lwd=0.5,coef=0,notch = 0.5,position = position_dodge(0.9)) +
  
  scale_fill_manual(values = palette[c(20,40)]) +
  scale_color_manual(values = palette[c(20,40)]) +
  
  # stat_summary(fun.y=mean, geom="point",aes(fill = warmt), 
  #              shape=23, size=2,show.legend = FALSE) +
  
  scale_x_discrete(labels=c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  scale_y_reverse(breaks = c(0,25,50,75,100),limits = c(100,0),labels=paste0(c(0,25,50,75,100),'%')) +
  xlab(label = ' ') +
  ylab(label = 'Range contraction') +
  
  theme_bw() +
  coord_cartesian(expand=F) +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(linetype = 'dashed',color='black'),
    axis.ticks.x = element_blank(),
    text = element_text(size=10),
    axis.text.x = element_text(color='black',vjust = 3),
    axis.text.y = element_text(color='black',angle=90,hjust = 0.5, vjust=1),
    axis.line.y = element_line(color='black'),
    # axis.line.y.right = element_line(),
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_text(angle = 0, size = 16)
    
    
  )

