source('R/4targets/MASTER.R');

library(raster); library(foreach); library(sf)
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

list.rasters.name <- paste0(dir_mod,dir_('figs/data/'),'raster_layers_figure_1.rds')

if(file.exists(list.rasters.name)){
  list.rasters <- readRDS(list.rasters.name)
}else{
  list.rasters <- foreach(n = seq_along(warming_targets)) %do% {
    
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
  names(list.rasters) <- warming_targets
  
  saveRDS(file = list.rasters.name,object = list.rasters)
  
}

list.rasters.name.dispersal <- paste0(dir_mod,dir_('figs/data/'),'raster_layers_dispersal_figure_1.rds')

if(file.exists(list.rasters.name)){
  list.rasters.dispersal <- readRDS(list.rasters.name.dispersal)
}else{
  list.rasters.dispersal <- foreach(n = seq_along(warming_targets)) %do% {
    
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
  names(list.rasters.dispersal) <- warming_targets
  
  saveRDS(file = list.rasters.name.dispersal,object = list.rasters.dispersal)
  
}

# start ~ 2:20PM

#> MAPVIEW ----------------------------------------------------------------------------------------------------
# library(mapview); library(rgdal); library(viridis)
# # add initial occurrence data
# SR <- readRDS(paste0(dir_mod,'gfdl/SR_tab_rcp4p5_1.5C_2044.rds'))
# SR <- SR[!is.na(SR@data$occ),]
# sr <- rasterize(SR,raster(ext=extent(SR),resolution=1/12),field = 'occ')
# # srp <- rasterToPolygons(sr,dissolve = T)
# # mapview(srp)
# # writeOGR(srp,paste0(dir_mod,'figs/data/sr_poly.gpkg'),'sr',driver = 'GPKG')
# 
# r <- brick(list.rasters)
# #try
# brksSR <- c(1,5,10,20,30,50,100,150,200,250,350)
# brksUniv <- c(seq(0,1,0.025))
# pal <- readRDS(paste0('pal1_40.rds')) #function
# #,maxpixels = 7212730
# mv <- mapview(sr,col.regions = rev(viridis_pal()(length(brksSR) - 1)),at = brksSR,legend = T,maxpixels = 7212730,na.color = 'transparent')
# mapshot(mv,url = 'D:/fishsuit/fishsuit_completeRun_warming_4targets/figs/data/SR.html',
#         remove_controls = NULL, selfcontained = FALSE)
# mv <- mapview(r,col.regions = rev(pal(40)),at = brksUniv,legend = T,maxpixels = 7212730,na.color = 'transparent')
# mapshot(mv,url = 'D:/fishsuit/fishsuit_completeRun_warming_4targets/figs/data/rel_loss_warming_targets4.html',
#         remove_controls = NULL, selfcontained = FALSE)

#> MAPS -------------------------------------------------------------------------------------------------------

r <- list.rasters[[3]]
rbase <- raster('ldd.asc')
rbase[!is.na(rbase)] <- -1

res(rbase) <- res(r)
rbase <- crop(rbase,extent(r))
extent(rbase) <- extent(r)

r[is.na(r) & rbase==-1] <- -1

# first try with viridis palette type 'c' > plasma
# use ggplot instead of levelplot?
# try to convert raster to polygons?
# library(colorspace);
# pal <- choose_palette()
# saveRDS(pal,'pal_YR_40.rds')

library(rasterVis); library(viridis); library(cowplot);library(grid);library(RColorBrewer)

palette <- c(
  brewer.pal(9, "Greys")[3],
  colorRampPalette(brewer.pal(9, "YlOrRd"))(40)
)

myTheme <- BTCTheme()
myTheme$panel.background$col = brewer.pal('RdBu',n=11)[11] #white

brksUniv <- c(-1,-0.1,seq(0.025,1,0.025))
ckey <- list(space='bottom',
             cex = 1.3
             ,col = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)
             ,at = seq(0,1,0.025)
             # ,labels =seq(0,1,0.2)
             ) #,labels=list(cex=1), height=1 


lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
)

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
            # ,margin = T
  )
)
p_main <- grid.grab()

dev.off()

for(n in c(1,2,4)){
  
  r <- list.rasters[[n]]
  r[is.na(r) & rbase==-1] <- -1
  
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
              ,colorkey=NULL
              
    )
  )
  p <- grid.grab()
  assign(paste0('p_panels',n),p)
  dev.off()
}

#> VIOLIN PLOTS INSET -----------------------------------------------------------------------------------------

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
        ESH_median = median(t$ESH,na.rm=T),
        warmt = wt
      )
    )
  } 
}

library(ggplot2)
p_inset <- ggplot(tab_median,aes(x=warmt,y=ESH_median)) + #
  geom_violin(aes(fill = warmt),draw_quantiles = c(0.25,0.5,0.75),lwd = .5,color='white') + #
  scale_fill_manual(values = palette[c(10,20,30,40)]) +
  # geom_boxplot(width = 0.08,outlier.size = 0.3,outlier.alpha = 0.5) +
  # geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  scale_y_reverse(breaks = c(0,50,100),limits = c(100,0),labels=paste0(c(0,50,100),'%')) +
  xlab(label = ' ') +
  # ylab(label = '% of suitable habitat lost') +
  ylab(' ') +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(linetype = 'dashed',color='white'),
    axis.ticks.x = element_blank(),
    text = element_text(size=10),
    axis.text.x = element_text(color='white',vjust = 3),
    axis.text.y = element_text(color='white',angle=90,hjust = 0.5, vjust=1),
    axis.line.y = element_line(color='white'),
    # axis.line.y.right = element_line(),
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    
  ) 
p_inset

library(grid); library(ggplotify); library(cowplot); library(ggpubr)

p_main_w_inset <- as.ggplot(p_main) + annotation_custom(grob = ggplotGrob(p_inset),xmin = 0.04, xmax = .28, 
                                                        ymin = .09, ymax = .48)

brow <- ggarrange(p_panels1,p_panels2,p_panels4,ncol=3,nrow=1,labels=letters[2:4],align = 'hv')
fig1 <- ggarrange(ggplotGrob(p_main_w_inset),brow,labels = c('a',''),ncol=1,nrow=2,heights = c(2.3, 1))

cowplot::ggsave(paste0(dir_mod,'figs/Figure_1_habitat_losses_warming_targets_Winset_alt2.jpg'),fig1,width = 220,height = 140,units='mm',dpi = 600)


