source('R/4targets/MASTER.R')

# Library
library(fmsb); library(foreach); library(grid)

#------------------------------------------------------------
#> DATA

# this is available from 'Figure_SI_boxplot_scenario.R'
tab <- readRDS(paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds'))
tab <- droplevels( tab[tab$ESH_type == 'Total',] )

ncores = 4
# parallelized ---
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
tab_median <- foreach(sp = unique(tab$id_no),.combine = 'rbind',.packages = 'foreach') %dopar% {
  tsp <- droplevels(tab[tab$id_no == sp,])
  foreach(wt = levels(tab$warmt),.combine = 'rbind') %do% {
    
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
parallel::stopCluster(cl)
# ---

tab_median4 <- droplevels(tab_median)
# retrieve IUCN metadata
iucn <- foreach(i = 1:2,.combine = 'rbind') %do% foreign::read.dbf(paste0('data/FW_FISH_PART_',i,'.dbf'))
iucn_simplified <- droplevels(foreach(i = unique(tab_median4$id_no),.combine = 'rbind') %do% iucn[iucn$id_no == i,][1,])

# habitat metadata
hab_type <- read.csv('data/iucn_habitat_type.csv')
hab_data <- foreach(i = unique(hab_type$binomial),.combine = 'rbind') %do% hab_type[hab_type$binomial == i,][1,] #lookslike there are more rows than species

# fishbase metadata
fishbase <- read.csv('data/iucn_fishbase.csv')

#hydro-climatic spatial features metadata
hydro <- read.csv('data/species_climate_feow_waterhed_area.csv')



#---------------------------------------------
#> CATEGORIES ADJUSTMENT

# create table
tab <- merge(tab_median4,iucn_simplified,by='id_no')

# IUCN code---
tab$code[tab$code == 'LR/lc'] <- 'LC'
tab$code[tab$code == 'EX'] <- 'DD'
# levels(tab$code) <- c('Critically Endangered','Data Deficient','Endangered',
#                       'Extinct','Least Concern','Near Threatened','Vulnerable','LR/lc')

tab <- merge(tab,hab_data[,c('binomial','lotic','lentic','seasonal')],by='binomial',all.x=T) # we lose 76 species when merging <<<<<<

# IUCN habitat type---
tab$habtype <- 'Le'
tab$habtype[tab$lotic == 1 & tab$lentic == 0] <- 'Lo'
tab$habtype[tab$lotic == 1 & tab$lentic == 1] <- 'LL'
tab$habtype[tab$marine == 't'] <- 'M' # consider any fish that is also marine as diadromous


tab <- droplevels(merge(tab,hydro,by = 'id_no'))

tab$rangesizecat <- '<4' #100 cells
tab$rangesizecat[tab$area >= 10**4] <- '4-4.5'
tab$rangesizecat[tab$area >= 10**4.5] <- '4.5-5'
tab$rangesizecat[tab$area >= 10**5] <- '5-6'
tab$rangesizecat[tab$area >= 10**6] <- '>6'
tab$rangesizecat <- factor(tab$rangesizecat,levels = c('>6','5-6','4.5-5','4-4.5','<4'))


# abbreviate names for order_
nchar_names <- 5
min_sample_size <- 20
levels(tab$order_) <- paste0(as.character(sapply(levels(tab$order_),function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')
minor_orders <- names(table(tab$order_)[which(table(tab$order_) <= min_sample_size*4)])
levels(tab$order_)[levels(tab$order_) %in% minor_orders] <- 'OTHER'





# koppen-geiger climatic zoning---
# plot(raster('data/Koeppen-Geiger-Classification-Reclassfied_3min_major.nc'))
# correct for 24 species with same of cells on 2 zones
tab$kg_main[tab$kg_main == '1;2'] <- 1
tab$kg_main[tab$kg_main == '1;3'] <- 1
tab <- droplevels(tab)
# assign names to levels
levels(tab$kg_main) <- LETTERS[1:5]
tab$kg_main <- factor(tab$kg_main,levels = rev(LETTERS[1:5]))

# levels(tab$kg_main) <- c(
#   'Equatorial','Arid','Warm-temperate','Snow','Polar'
# )
#A: Equatorial
#B: Arid
#C: Warm-temperate
#D: Snow
#E: Polar

# # Ecoregions realms---
# tab$feow_main[tab$feow_main == '2;3'] <- 2
# tab <- droplevels(tab)
# levels(tab$feow_main) <- c(
#   'Nearctic','Neotropic','Neotropic','Palearctic','Afrotropic','Palearctic','Indo-Malay','Australasia','NA'
# )


tab <- droplevels(merge(tab,fishbase,by.x = 'binomial',by.y = 'iucn_name',all.x=T))


#Migratory---
#Anacat #<<<<<<<<< there are about 52 oceanodromous fishes (spend entire life in ocean, weird!!)
levels(tab$AnaCat) <- c(NA,rep('Diad.',5),'Non.',NA,'Pota.')

#Importance---
levels(tab$Importance) <- c(NA,'Com.','HCom.','MCom.','NoInt.','NoInt.','Subs.')

#PriceCateg---
levels(tab$PriceCateg) <- c('hi','lo','med',NA,'vhi')

#Body Length
tab$lengthcat[!is.na(tab$Length)] <- '<5'
tab$lengthcat[tab$Length >= 5] <- '5-10'
tab$lengthcat[tab$Length >= 10] <- '10-30'
tab$lengthcat[tab$Length >= 30] <- '30-80'
tab$lengthcat[tab$Length >= 80] <- '80-200'
tab$lengthcat[tab$Length >= 200] <- '>200'
tab$lengthcat <- factor(tab$lengthcat,levels = c('>200','80-200','30-80','10-30','5-10','<5'))


#FoodTroph
tab$foodtrophcat[!is.na(tab$FoodTroph)] <- 'Herbi.'
tab$foodtrophcat[tab$FoodTroph > 2.19 & tab$FoodTroph <= 2.79] <- 'Omni.'
tab$foodtrophcat[tab$FoodTroph > 2.79] <- 'Carni.'

#Vulnerability #interesting to compare with the IUCN categories, e.g. boxplots
#Anacat #train an ANN to validate migratory categories inferred from IUCN with fishbase data, maybe phylogenetic approach?

tab <- droplevels(tab)

#---------------------------------------------
#> RADAR PLOTS
library(ggplot2); library(grid); library(ggplotify); library(RColorBrewer)

catv <- c('kg_main','habtype','AnaCat','rangesizecat','lengthcat','foodtrophcat','code','Importance','order_')
# nchar_namesv <- c(rep(NA,5),4)
tit = NULL

# n_minv=c(rep(5,5),20)
size_namesv = c(rep(0.8,8),0.5)

radplots <- foreach(i = seq_along(catv)) %do% {
  
  cat <- catv[i]
  n_min <- 1 #n_minv[i]
  nchar_names <- NA #nchar_namesv[i]
  size_names <- size_namesv[i]
  
  data <- foreach(wt = levels(tab$warmt),.combine = 'rbind') %do% {
    t <- tab[tab$warmt == wt,]
    ts <- split(t,t[,cat])
    res <- do.call('cbind',lapply(ts,function(x) {x = x[!is.na(x$ESH_median),]; c(nrow(x),sum(x$ESH_median > 50)/nrow(x)*100 )}))
    row.names(res) <- c('n',wt)
    if(wt == '1.5C'){
      return(as.data.frame(res))
    }else{
      res2 <- t(as.data.frame(res[2,]))
      row.names(res2) <- wt
      return(res2)
    }
  }
  
  data <- data[,data[1,] >= n_min]
  
  new_names <- colnames(data)
  if(!is.na(nchar_names)) new_names <- paste0(as.character(sapply(colnames(data),function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')
  
  colnames(data) <- paste0(new_names,'\n(',as.integer(data['n',]),')')
  # if(!is.na(nchar_names)) colnames(data) <- paste0(new_names,'\n(',as.integer(data['n',]),')')
  
  data <- rbind(rep(100,ncol(data)),rep(0,ncol(data)),data[2:nrow(data),][rev(row.names(data[2:nrow(data),])),])
  
  # custom_pal <- pal(5)[2:5]
  custom_pal <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)])
  colors_border<- custom_pal
  colors_in <- do.call('c',lapply(custom_pal,function(x) rgb(t(col2rgb(x)) ,alpha = 180,maxColorValue = 255)))
  
  
  # par(mar = c(0,0,0,0),oma = c(0,0,0,0))
  p <- ggplotGrob(
    as.ggplot(
      ~radarchart( data  , axistype=1 , pty=32, seg=4,
                   #custom polygon
                   pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
                   #custom the grid
                   cglty=1, cglwd=1, cglcol="black", axislabcol="black", calcex = 0.7, caxislabels=c('0%','','50%','','100%'), #caxislabels=paste0(seq(0,100,25),'%'),
                   #custom labels
                   vlcex=size_names,title = tit)
    ) + theme(plot.margin = unit(c(-0.5, -0.5, -2, -2), "cm"))
  )
  # if(legendp) print(legend(x=1.1, y=1.4, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3))
  
  return(p)
  
}

# plot.new()
# leg <- as.grob(
#  
#   ~legend(x=1.1, y=1.4, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
# )
#---------------------------------------------
#> SCATTER PLOTS

# alphas <- 0.3 #dots alpha
# fitmethod <- 'auto' #'lm'
# 
# # continuous plots
# library(ggplot2); library(RColorBrewer)
# p_area <- ggplot(tab,aes(x = log10(area),y = ESH_median,color = warmt,fill = warmt)) +
#   geom_point(alpha = alphas) +
#   geom_smooth(method = fitmethod,color=brewer.pal(9, "Greys")[7],alpha = 0.5) +
#   scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)]) +
#   scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)]) +
#   xlab('LOG-Range area (km2)') +
#   ylab('Range lost (%)') +
#   coord_cartesian(ylim = c(0,101),expand = FALSE) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),legend.position = "none",
#         panel.border = element_blank(),
#         axis.line = element_line())
# 
# p_length <- ggplot(tab,aes(x = log10(Length),y = ESH_median,color = warmt,fill = warmt)) +
#   geom_point(alpha = alphas) +
#   geom_smooth(method = fitmethod,color=brewer.pal(9, "Greys")[7],alpha = 0.5) +
#   scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)]) +
#   scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)]) +
#   xlab('LOG-Body length (cm)') +
#   ylab(' ') +
#   coord_cartesian(ylim = c(0,101),expand = FALSE) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),legend.position = "none",
#         panel.border = element_blank(),
#         axis.line = element_line())
# 
# p_troph <- ggplot(tab,aes(x = FoodTroph,y = ESH_median,color = warmt,fill = warmt)) +
#   geom_point(alpha = alphas) +
#   geom_smooth(method = fitmethod,color=brewer.pal(9, "Greys")[7],alpha = 0.5) +
#   scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)]) +
#   scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)]) +
#   xlab('Trophic level (-)') +
#   ylab(' ') +
#   coord_cartesian(ylim = c(0,101),expand = FALSE) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),legend.position = "none",
#         panel.border = element_blank(),
#         axis.line = element_line())
# 
# 
# scatterplots <- list(
#   ggplotGrob(p_area),ggplotGrob(p_length),ggplotGrob(p_troph)
# )

#------------------------------------------------------------
#> LEGEND
p_leg <- ggplot(tab,aes(x = log10(area),y = ESH_median,color = warmt)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                     labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  guides(colour = guide_legend(title=NULL,override.aes = list(size = 3))) +
  theme_bw() +
  theme(legend.direction = 'horizontal',
        text = element_text(size = 15))

tmp <- ggplot_gtable(ggplot_build(p_leg)) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 


# legend <- g_legend(p_leg) 
grid.draw(legend)


#------------------------------------------------------------
#> DRAWING
# library(ggpubr)
# fig <- ggarrange(
#   ggarrange(plotlist=radplots,ncol= 3, nrow = 2, labels = LETTERS[1:6]),
#   legend,
#   ggarrange(plotlist=scatterplots,ncol= 3, nrow = 1, labels = LETTERS[7:9]),
#   ncol = 1,nrow = 3,heights = c(2,0.1,1))
# cowplot::ggsave(paste0(dir_mod,'figs/Figure_4_maj50_1.jpg'),fig,width = 220,height = 230,units='mm',dpi = 600)

library(ggpubr)
fig <- ggarrange(
  ggarrange(plotlist=radplots,ncol= 3, nrow = 3, labels = letters[1:9]),
  legend,
  ncol = 1,nrow = 2,heights = c(3,0.1))
cowplot::ggsave(paste0(dir_mod,'figs/Figure_3.jpg'),fig,width = 220,height = 230,units='mm',dpi = 600)

