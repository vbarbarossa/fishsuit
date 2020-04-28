source('R/4targets/MASTER.R');

library(raster)

tab <- readRDS(paste0(dir_mod,'hadgem/SR_tab_rcp8p5_3.2C_2058.rds'))
tab <- tab[!is.na(tab@data$occ),]

ras <- rasterize(tab,raster(ext=extent(tab),resolution=1/12),field = 'occ') #can take a few minutes

#--------------------------------------------------------------------------------------
#>> main map
library(rasterVis); library(viridis); library(RColorBrewer); library(cowplot);library(grid)

# raster with number of species per grid-cell
r <- ras

rbase <- raster('ldd.asc')
rbase[!is.na(rbase)] <- -1

res(rbase) <- res(r)
rbase <- crop(rbase,extent(r))
extent(rbase) <- extent(r)

r[is.na(r) & rbase==-1] <- -1


brksUniv <- c(-1,0,5,10,20,30,50,100,150,200,250,350)

palette <- c(
  brewer.pal(9, "Greys")[3],
  rev(viridis_pal()(length(brksUniv) - 2))
)

myTheme <- BTCTheme()
myTheme$panel.background$col = 'white'

ckey <- list(space='bottom',
             cex = 1.3
             ,col = rev(viridis_pal()(length(brksUniv) - 2))
             ,at = c(1,5,10,20,30,50,100,150,200,250,350),
             labels=list(at = c(1,10,30,50,100,150,200,250,350),
                         labels=c(1,10,30,50,100,150,200,250,350))
)

# ckey <- #list(space='bottom',cex = 1.3) #,labels=list(cex=1), height=1
#   list(
#     space = 'bottom',
#     cex = 1.3,
#     labels=list(at = brksUniv, 
#                 labels=brksUniv
#     )
#   )

lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=1.5), right.padding=list(x=1.5))
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

p <- grid.grab()

dev.off()

cowplot::ggsave(paste0(dir_mod,'figs/Figure_SI_original_SR.jpg'),p,width = 220,height = 120,units='mm',dpi = 600)


#------------------------------------------------------------------------------------------------------------------
# lotic vs lentic

tab <- readRDS(paste0(dir_mod,'hadgem/SR_tab_rcp8p5_3.2C_2058_lotic.rds'))
tab <- tab[!is.na(tab@data$occ),]

ras_lo <- rasterize(tab,raster(ext=extent(tab),resolution=1/12),field = 'occ') #can take a few minutes

tab <- readRDS(paste0(dir_mod,'hadgem/SR_tab_rcp8p5_3.2C_2058_lentic.rds'))
tab <- tab[!is.na(tab@data$occ),]

ras_le <- rasterize(tab,raster(ext=extent(tab),resolution=1/12),field = 'occ') #can take a few minutes


#--------------------------------------------------------------------------------------
plot_sr <- function(r,legend = TRUE){
  
  rbase <- raster('ldd.asc')
  rbase[!is.na(rbase)] <- -1
  
  res(rbase) <- res(r)
  rbase <- crop(rbase,extent(r))
  extent(rbase) <- extent(r)
  
  r[is.na(r) & rbase==-1] <- -1
  
  
  brksUniv <- c(-1,0,5,10,20,30,50,100,150,200,250,350)
  
  palette <- c(
    brewer.pal(9, "Greys")[3],
    rev(viridis_pal()(length(brksUniv) - 2))
  )
  
  myTheme <- BTCTheme()
  myTheme$panel.background$col = 'white'
  
  ckey <- list(space='bottom',
               cex = 1.3
               ,col = rev(viridis_pal()(length(brksUniv) - 2))
               ,at = c(1,5,10,20,30,50,100,150,200,250,350),
               labels=list(at = c(1,10,30,50,100,150,200,250,350),
                           labels=c(1,10,30,50,100,150,200,250,350))
  )
  
  lattice.options(
    layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
    layout.widths=list(left.padding=list(x=1.5), right.padding=list(x=1.5))
  )
  
  if(!legend) ckey <- NULL
  
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
  
  p <- grid.grab()
  
  return(p)
  
}

p_lo <- plot_sr(ras_lo,legend = FALSE)
p_le <- plot_sr(ras_le,legend = TRUE)

cowplot::ggsave(paste0(dir_mod,'figs/Figure_SI_original_SR_lentic_vs_lotic.jpg'),
                ggpubr::ggarrange(p_lo,p_le,labels = c('Lotic','Exclusively lentic'),hjust = c(-0.65,-0.2),ncol=1,nrow=2,heights = c(1,1.1),align = 'hv'),
                width = 160,height = 150,units='mm',dpi = 600)

