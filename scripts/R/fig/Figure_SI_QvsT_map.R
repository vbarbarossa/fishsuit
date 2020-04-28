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
rasterize_rel_losses <- function(x,var,dir_root = 'fishsuit_completeRun_warming_4targets'){
  tab <- readRDS(paste0(dir_root,'/',atab$clmod[x],'/SR_tab_',
                        atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
  tab <- tab[!is.na(tab@data$occ),]
  tab@data <- data.frame(rel = (tab@data$occ - tab@data[,var])/tab@data$occ)
  return(
    rasterize(tab,raster(ext=extent(tab),resolution=1/12),field = 'rel')
  )
}


varnames <- c('Q_all','T_all','both_QT')

n = 3 #3.2C scen
list.rasters.name <- paste0(dir_mod,dir_('figs/data/'),'raster_layers_figure_2_',warming_targets[n],'.rds')
if(file.exists(list.rasters.name)){
  list_3C <- readRDS(list.rasters.name)
}else{
  list_3C <- foreach(vn = varnames) %do% {
    
    atab <- attr.tab[attr.tab$warmt == warming_targets[n],]
    
    # rasterize_rel_losses(x=1,var = vn)
    
    ncores = 7
    # parallelized ---
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    br <- brick(
      foreach(i = 1:nrow(atab),.packages='raster') %dopar% rasterize_rel_losses(x=i,var = vn)
    )
    parallel::stopCluster(cl)
    # ---
    return(
      calc(br,function(x) median(x,var=vn,na.rm=T))
    )
    
  }
  names(list.rasters) <- varnames
  
  saveRDS(file = list.rasters.name,object = list_3C)
  
}
#> MAPS -------------------------------------------------------------------------------------------------------

# choose palette for different variables
# library(colorspace);
# pal <- choose_palette()
# saveRDS(pal,'pal_QT.rds')

# palnames <- c('Q','T','QT')
# 
# library(rasterVis); library(viridis); library(cowplot);library(grid)
# 
# myTheme <- BTCTheme()
# myTheme$panel.background$col = 'white'
# 
# ckey <- list(space='bottom',cex = 1.3) #,labels=list(cex=1), height=1 
# brksUniv <- c(seq(0,1,0.025))
# 
# lattice.options(
#   layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
#   layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
# )

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


r <- list_3C[[1]]
rbase <- raster('ldd.asc')
rbase[!is.na(rbase)] <- -1
res(rbase) <- res(r)
rbase <- crop(rbase,extent(r))
extent(rbase) <- extent(r)

pgrobs <- foreach(n = 1:3) %do% {
  
  r <- list_3C[[n]]
  r[is.na(r) & rbase==-1] <- -1
  
  # pal <- readRDS(paste0('pal1_40.rds')) #function
  # # pal <- readRDS(paste0('pal_',palnames[n],'.rds')) #function
  # palette <- rev(pal(40))
  
  if(n %in% 1:2){
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
                # ,margin = T
      )
    )
    p <- grid.grab()
    
  }else{
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
                ,colorkey=ckey #this changes to ckey if want to display legend on last panel only
                # ,margin = T
      )
    )
    p <- grid.grab()
    
  }
  
  dev.off()
  
  return(p)
  
}



tab <- readRDS(paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds'))
tab <- droplevels(tab[tab$warmt == paste0(warming_targets[n],'C'),])

tab_median_vars <- foreach(sp = unique(tab$id_no),.combine = 'rbind') %do% {
  tsp <- droplevels(tab[tab$id_no == sp,])
  foreach(vn = levels(tab$ESH_type),.combine = 'rbind') %do% {
    
    t <- droplevels(tsp[tsp$ESH_type == vn,])
    return(
      data.frame(
        id_no = sp,
        no.cells = unique(t$no.cells),
        ESH_median = median(t$ESH,na.rm=T),
        varname = vn
      )
    )
  } 
}

tab_median_vars$ESH_weighted <- tab_median_vars$no.cells*tab_median_vars$ESH_median

T_loss <- round(sum(tab_median_vars$ESH_weighted[tab_median_vars$varname == 'Tw'],na.rm=T)*100/10**9,1) #G km2
Q_loss <- round(sum(tab_median_vars$ESH_weighted[tab_median_vars$varname == 'Q'],na.rm=T)*100/10**9,1) #G km2
QT_loss <- round(sum(tab_median_vars$ESH_weighted[tab_median_vars$varname == 'Q&Tw'],na.rm=T)*100/10**9,1) #G km2

library(eulerr); library(RColorBrewer) #<<<<<<<<<<< look into VennDiagram package (seems to have nicer aestetics)
con <- c(A = T_loss, B = Q_loss, "A&B" = QT_loss)
print(
  plot(euler(con), labels = c('Tw','Q'),
       fill = RColorBrewer::brewer.pal(n = 9,'Greys')[c(2,4)],
       # fill = c(rgb(228,26,28,alpha = 200,maxColorValue = 255),
       # rgb(55,126,184,alpha = 200,maxColorValue = 255)),
       col = rep('white',2),
       quantities = T)
)
peuler <- grid.grab()

fig <- ggpubr::ggarrange(pgrobs[[1]],pgrobs[[2]],pgrobs[[3]],peuler,labels = letters[1:4],ncol=1,nrow=4,heights = c(1,1,1.15,0.5))

cowplot::ggsave(paste0(dir_mod,'figs/Figure_SI_QT_contribution.jpg'),fig,width = 220,height = 350,units='mm',dpi = 600)



# palQ <- rev(readRDS('pal_Q.rds')(5))[4]
# palT <- rev(readRDS('pal_T.rds')(5))[4]
# palQT <- rev(readRDS('pal_QT.rds')(5))[4]
# pal_all <- rev(readRDS('pal1_40.rds')(5))[4]

# library(ggplot2)
# p <- ggplot(tab_median_vars,aes(x=varname,y=ESH_median)) + #
#   geom_violin(aes(color = varname),draw_quantiles = T,lwd = .5, scale = 'width',fill = 'transparent') + #
#   geom_boxplot(width = 0.1,outlier.size = 1,outlier.alpha = 0.5) +
#   scale_color_manual(values = c(pal_all,palQ,palT,palQT)) +
#   scale_x_discrete(labels=c('Total','Flow','Water temp.','Overlap')) +
#   scale_y_reverse(breaks = c(0,50,100),limits = c(100,0),labels=paste0(c(0,50,100),'%')) +
#   xlab(label = ' ') +
#   # ylab(label = '% of suitable habitat lost') +
#   ylab(' ') +
#   theme_bw() +
#   theme(
#     legend.position="none",
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     panel.grid.major.y = element_line(linetype = 'dashed'),
#     axis.ticks.x = element_blank(),
#     text = element_text(size=12),
#     axis.text.x = element_text(color='black',vjust = 3),
#     axis.text.y = element_text(color='black',angle=90,hjust = 0.5, vjust=1),
#     axis.line.y = element_line(),
#     # axis.line.y.right = element_line(),
#     panel.background = element_rect(fill = "transparent") # bg of the panel
#     , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
#     , panel.grid.major = element_blank() # get rid of major grid
#     , panel.grid.minor = element_blank() # get rid of minor grid
#     , legend.background = element_rect(fill = "transparent") # get rid of legend bg
#     , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
#     
#   ) 
# p

p <- ggplot(droplevels(tab_median_vars[tab_median_vars$varname != 'Total',]),aes(x=varname,y=ESH_median)) + #
  # geom_violin(aes(fill = varname),draw_quantiles = T,lwd = .5, scale = 'area') + #
  geom_boxplot(width = 0.3,outlier.size = 1,outlier.alpha = 0.5) +
  # scale_fill_manual(values = c(palQ,palT,palQT)) +
  scale_x_discrete(labels=c('Flow','Water temp.','Overlap')) +
  scale_y_reverse(breaks = c(0,50,100),limits = c(100,0),labels=paste0(c(0,50,100),'%')) +
  xlab(label = ' ') +
  # ylab(label = '% of suitable habitat lost') +
  ylab(' ') +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(linetype = 'dashed'),
    axis.ticks.x = element_blank(),
    text = element_text(size=12),
    axis.text.x = element_text(color='black',vjust = 3),
    axis.text.y = element_text(color='black',angle=90,hjust = 0.5, vjust=1),
    axis.line.y = element_line(),
    # axis.line.y.right = element_line(),
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    
  ) 
# p

fig <- ggpubr::ggarrange(pgrobs[[1]],pgrobs[[2]],pgrobs[[3]],ggplotGrob(p),labels = letters[1:4],ncol=2,nrow=2)

cowplot::ggsave(paste0(dir_mod,'figs/Figure_2_QT_contribution_alt2.jpg'),fig,width = 220,height = 110,units='mm',dpi = 600)



# -------------------------
# do euler plots at each warming target

T_loss <- numeric()
Q_loss <- numeric()
QT_loss <- numeric()

for(n in 1:4){
  
  tab <- readRDS(paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds'))
  tab <- droplevels(tab[tab$warmt == paste0(warming_targets[n],'C'),])
  
  tab_median_vars <- foreach(sp = unique(tab$id_no),.combine = 'rbind') %do% {
    tsp <- droplevels(tab[tab$id_no == sp,])
    foreach(vn = levels(tab$ESH_type),.combine = 'rbind') %do% {
      
      t <- droplevels(tsp[tsp$ESH_type == vn,])
      return(
        data.frame(
          id_no = sp,
          no.cells = unique(t$no.cells),
          ESH_median = median(t$ESH,na.rm=T),
          varname = vn
        )
      )
    } 
  }
  
  tab_median_vars$ESH_weighted <- tab_median_vars$no.cells*tab_median_vars$ESH_median
  
  T_loss[n] <- round(sum(tab_median_vars$ESH_weighted[tab_median_vars$varname == 'Tw'],na.rm=T)*100/10**9,1) #G km2
  Q_loss[n] <- round(sum(tab_median_vars$ESH_weighted[tab_median_vars$varname == 'Q'],na.rm=T)*100/10**9,1) #G km2
  QT_loss[n] <- round(sum(tab_median_vars$ESH_weighted[tab_median_vars$varname == 'Q&Tw'],na.rm=T)*100/10**9,1) #G km2
  
}

library(eulerr); library(RColorBrewer)
con <- c(A = T_loss[1], B = Q_loss[1], "A&B" = QT_loss[1],
         C = T_loss[2], D = Q_loss[2], "C&D" = QT_loss[2],
         E = T_loss[3], F = Q_loss[3], "E&F" = QT_loss[3],
         G = T_loss[4], H = Q_loss[4], "G&H" = QT_loss[4])

jpeg(paste0(dir_mod,'figs/Figure_2_QT_contribution_eulerONLY.jpg'),width = 220,height = 110,units='mm',res = 600)
plot(euler(con), 
     labels = NULL,
     fills = rep(c(rgb(228,26,28,alpha = 200,maxColorValue = 255),
                   rgb(55,126,184,alpha = 200,maxColorValue = 255)),4),
     quantities = T)

dev.off()
