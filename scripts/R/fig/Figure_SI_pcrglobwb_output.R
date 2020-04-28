source('R/4targets/MASTER.R');

source_layers <- paste0(dir_mod,'figs/data/pcrglobwb_anomaly_ensemble.rds')

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
    foreach(var = vars) %do% paste0(dir_mod,climate_models,'/pcrglobwb_processed/merged/',var,'_hist.tif'),
    ensemble_median
  )
  names(hist) <- vars
  
  diff_ens <- foreach(wt = unique(attr.tab$warmt)) %do% {
    
    attrs <- droplevels(attr.tab[attr.tab$warmt == wt,])
    
    fut <- lapply(
      foreach(var = vars) %do% {
        foreach(r = 1:nrow(attrs),.combine = 'c') %do% paste0(dir_mod,attrs[r,'clmod'],
                                                              '/pcrglobwb_processed/merged/',var,'_',
                                                              attrs[r,'scen'],'_',attrs[r,'warmt'],'C_',
                                                              attrs[r,'year'],'.tif')
      },
      ensemble_median
    )
    names(fut) <- vars
    
    diff <- lapply(1:3,function(x) {mask(fut[[x]],hist[[x]] >= 0) - mask(hist[[x]],hist[[x]] >= 0)})
    names(diff) <- vars
    
    return(diff)
    
  }
  
  names(diff_ens) <- warming_targets
  saveRDS(diff_ens,paste0(dir_mod,'figs/data/pcrglobwb_anomaly_ensemble.rds'))
  
}

#--------------------------------------------------------------------------------------
# one plot for each variable at 4 different warming targets

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
df <- rbind(
  foreach(i = seq_along(deg),.combine='rbind') %do% {
    as(diff_ens[[warming_targets[i]]][['Qzf']] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(warming = deg[i]) %>%
      dplyr::select(x,y,value = layer,warming)
  }
)

df$warming <- factor(df$warming, levels = deg)

# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_tile(data = df, aes(x=x, y=y, fill=value), alpha=0.8) +
  # scale_fill_viridis_c(breaks = seq(0,1,0.1),
  #                      labels = seq(0,1,0.1),
  #                      limits = c(0,1),
  #                      option = 'C',na.value = "transparent") +
  facet_grid(warming~.) +
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
ggsave(paste0(dir_mod,'figs/maps_pcrglobwb_anomaly_Qzf.jpg'),p,
       width = 183,height = 200,dpi = 600,units = 'mm')






library(rasterVis); library(viridis); library(RColorBrewer); library(cowplot);library(grid)
library(colorspace); 
# pal <- choose_palette()
# saveRDS(pal,'palSIpcrglobwb_output.rds')
pal <- readRDS('palSIpcrglobwb_output.rds') #function

myTheme <- BTCTheme()
myTheme$panel.background$col = 'white'

lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=1.5), right.padding=list(x=1.5))
)

ck <- list(space='bottom',cex = 1.3) #,labels=list(cex=1), height=1 

rbase <- raster('ldd.asc')
rbase[!is.na(rbase)] <- 1

#Qmi
dl <- lapply(diff_ens,function(x) x[[1]])
brksUniv <- c(-8000,-5000,-3000,-1000,-500,-100,-50,-20,-10,-5,-0.5,0.5,5,10,20,50,100,500,1000,3000,5000,8000)
palette <- rev(pal((length(brksUniv) - 1)))

r <- dl[[1]]
rbase <- raster('ldd.asc')
rbase[!is.na(rbase)] <- 1
res(rbase) <- res(r)
rbase <- crop(rbase,extent(r))
extent(rbase) <- extent(r)


l <- foreach(i = 1:4) %do% {
  r <- dl[[i]]
  
  r[is.na(rbase)] <- NA
  
  ckey <- NULL
  if(i == 4) ckey <- ck
  
  print(
    levelplot(r,
              maxpixels=1e10,
              margin=list(draw=FALSE),
              par.settings = myTheme,
              col.regions = palette,
              xlab=NULL,
              ylab=NULL,
              scales=list(draw=FALSE),
              pretty=F
              ,at=brksUniv
              ,colorkey=ckey
    )
  )
  p <- grid.grab()
  dev.off()
  return(p)
  
}

cowplot::ggsave(paste0(dir_mod,'figs/Figure_SI_pcrglobwb_anomaly_Qmi.jpg'),
                ggpubr::ggarrange(plotlist=l,labels = paste0(warming_targets,'C'),ncol=1,nrow=4,heights = c(1,1,1,1.21),align = 'hv'),
                width = 160,height = 230,units='mm',dpi = 600)


#Qzf
dl <- lapply(diff_ens,function(x) x[[2]])
brksUniv <- c(-52,-40,-30,-20,-10,-5,-2,-0.5,0.5,2,5,10,20,30,40,52)
palette <- pal((length(brksUniv) - 1))

r <- dl[[1]]
rbase <- raster('ldd.asc')
rbase[!is.na(rbase)] <- 1
res(rbase) <- res(r)
rbase <- crop(rbase,extent(r))
extent(rbase) <- extent(r)


l <- foreach(i = 1:4) %do% {
  r <- dl[[i]]
  
  r[is.na(rbase)] <- NA
  
  ckey <- NULL
  if(i == 4) ckey <- ck
  
  print(
    levelplot(r,
              maxpixels=1e10,
              margin=list(draw=FALSE),
              par.settings = myTheme,
              col.regions = palette,
              xlab=NULL,
              ylab=NULL,
              scales=list(draw=FALSE),
              pretty=F
              ,at=brksUniv
              ,colorkey=ckey
    )
  )
  p <- grid.grab()
  dev.off()
  return(p)
  
}

cowplot::ggsave(paste0(dir_mod,'figs/Figure_SI_pcrglobwb_anomaly_Qzf.jpg'),
                ggpubr::ggarrange(plotlist=l,labels = paste0(warming_targets,'C'),ncol=1,nrow=4,heights = c(1,1,1,1.21),align = 'hv'),
                width = 160,height = 230,units='mm',dpi = 600)


#Tma
dl <- lapply(diff_ens,function(x) x[[3]])
brksUniv <- c(seq(-30,-5,5),-4,-3,-2,-1,-0.5,0.5,1,2,3,4,seq(5,30,5))
palette <- pal((length(brksUniv) - 1))

r <- dl[[1]]
rbase <- raster('ldd.asc')
rbase[!is.na(rbase)] <- 1
res(rbase) <- res(r)
rbase <- crop(rbase,extent(r))
extent(rbase) <- extent(r)


l <- foreach(i = 1:4) %do% {
  r <- dl[[i]]
  
  r[is.na(rbase)] <- NA
  
  ckey <- NULL
  if(i == 4) ckey <- ck
  
  print(
    levelplot(r,
              maxpixels=1e10,
              margin=list(draw=FALSE),
              par.settings = myTheme,
              col.regions = palette,
              xlab=NULL,
              ylab=NULL,
              scales=list(draw=FALSE),
              pretty=F
              ,at=brksUniv
              ,colorkey=ckey
    )
  )
  p <- grid.grab()
  dev.off()
  return(p)
  
}

cowplot::ggsave(paste0(dir_mod,'figs/Figure_SI_pcrglobwb_anomaly_Tma.jpg'),
                ggpubr::ggarrange(plotlist=l,labels = paste0(warming_targets,'C'),ncol=1,nrow=4,heights = c(1,1,1,1.21),align = 'hv'),
                width = 160,height = 230,units='mm',dpi = 600)

