source('R/4targets/MASTER.R')

# Library
library(fmsb); library(foreach); library(grid); library(RColorBrewer)

#------------------------------------------------------------
#> DATA

ESH_filename <- 'tabs/ESH_merged_warmingtargets_confidencemargin.rds'
if(file.exists(ESH_filename)){
  tab_long <- readRDS(ESH_filename)
}else{
  
  tab_long <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    t0 <- read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab.csv'))
    t <- cbind(t0,as.data.frame(do.call('rbind',strsplit(as.character(t0[,'comboscen']),'_'))))
    t$GCM <- clmod
    t$comboscen <- paste0(t$comboscen,'_',t$GCM)
    colnames(t)[8:10] <- c('RCP','warmt','year')
    
    return(t)
  }
  tab_long <- droplevels(tab_long)
  tab_long$GCM <- factor(tab_long$GCM)
  tab_long$comboscen <- factor(tab_long$comboscen)
  saveRDS(tab_long,ESH_filename)
}

library(tidyr)
tab <- spread(tab_long[,c(1,3,7)], key = comboscen, value = ESH_all)

# function that selects the columns based on warming target
sel_col <- function(t,wt) sapply(strsplit(colnames(t)[1:ncol(t)],'_'),function(x) paste0(wt,'C') %in% x)
#str(tab[,sel_col(tab,'1.5')])

# retrieve IUCN metadata
iucn <- foreach(i = 1:2,.combine = 'rbind') %do% foreign::read.dbf(paste0('data/FW_FISH_PART_',i,'.dbf'))
iucn_simplified <- droplevels(foreach(i = unique(tab$id_no),.combine = 'rbind') %do% iucn[iucn$id_no == i,][1,])

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
tab <- merge(tab,iucn_simplified,by='id_no')

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


# save complete list of names for SI
# record level now
# t4n <- tab
# t4n$order_orig <- t4n$order_
# nchar_names <- 5
# min_sample_size <- 20
# # levels(t4n$order_) <- paste0(as.character(sapply(levels(t4n$order_),function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')
# minor_orders <- names(table(t4n$order_)[which(table(t4n$order_) <= min_sample_size)])
# levels(t4n$order_)[levels(t4n$order_) %in% minor_orders] <- 'OTHER'
# which(!levels(t4n$order_orig) %in% levels(t4n$order_))
# dn <- data.frame(name = levels(t4n$order_orig))
# dn$abbr <- as.character(dn$name)
# dn$abbr <- paste0(as.character(sapply(dn$abbr,function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')
# dn$abbr[which(!levels(t4n$order_orig) %in% levels(t4n$order_))] <- 'OTHER'
# write.csv(dn,paste0(dir_mod,'tabs/Figure_3_list_of_order_names.csv'),row.names = F)

# abbreviate names for order_
nchar_names <- 7
min_sample_size <- 20
levels(tab$order_)[levels(tab$order_) %in% c("CYPRINIFORMES","CYPRINODONTIFORMES")] <- paste0(as.character(sapply(levels(tab$order_)[levels(tab$order_) %in% c("CYPRINIFORMES","CYPRINODONTIFORMES")],function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')
nchar_names <- 5
levels(tab$order_)[!levels(tab$order_) %in% c("CYPRINI.","CYPRINO.")] <- paste0(as.character(sapply(levels(tab$order_)[!levels(tab$order_) %in% c("CYPRINI.","CYPRINO.")],function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')
minor_orders <- names(table(tab$order_)[which(table(tab$order_) <= min_sample_size)])
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

catv <- c('kg_main','rangesizecat','lengthcat','habtype','AnaCat','code','foodtrophcat','Importance','order_')
# nchar_namesv <- c(rep(NA,5),4)
tit = NULL

# n_minv=c(rep(5,5),20)
size_namesv = c(rep(0.8,8),0.5)

radplots <- foreach(i = seq_along(catv)) %do% {
  
  cat <- catv[i]
  n_min <- 1 #n_minv[i]
  nchar_names <- NA #nchar_namesv[i]
  size_names <- size_namesv[i]
  
  data <- foreach(wt = warming_targets,.combine = 'rbind') %do% {
    #split based on categories
    ts <- split(tab,tab[,cat])
    # calculate the mean proportion of species losing >50% of habitat at each scenario
    res <- do.call('cbind',
                   lapply(ts,function(x) {
                     x = x[,sel_col(x,wt=wt)]
                     x = x[complete.cases(x),]
                     # calculate the proportion of species losing >50% of habitat for each scenario
                     # then average the values across the scenarios
                     c(nrow(x),mean(apply(x,2,function(c) {sum(c > 50)/length(c)*100})) )
                   }))
    row.names(res) <- c('n',wt)
    if(wt == '1.5'){
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
  
  data <- rbind(rep(100,ncol(data)),rep(0,ncol(data)),data[2:nrow(data),][rev(row.names(data[2:nrow(data),])),])
  
  custom_pal <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)])
  colors_border<- custom_pal
  colors_in <- do.call('c',lapply(custom_pal,function(x) rgb(t(col2rgb(x)) ,alpha = 180,maxColorValue = 255)))
  
  
  # par(mar = c(0,0,0,0),oma = c(0,0,0,0))
  library(ggplot2); library(ggplotify)
  p <- ggplotGrob(
    as.ggplot(
      ~radarchart( data  , axistype=1 , pty=32, seg=4,
                   #custom polygon
                   pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
                   #custom the grid
                   cglty=1, cglwd=1, cglcol="black", axislabcol="black", calcex = 0.7, caxislabels=c('0%','','50%','',''), #caxislabels=paste0(seq(0,100,25),'%'),
                   #custom labels
                   vlcex=size_names,title = tit)
    ) + theme(plot.margin = unit(c(-0.5, -0.5, -2, -2), "cm"))
  )
  
  return(p)
  
}

#------------------------------------------------------------
#> LEGEND
p_leg <- ggplot(tab_long,aes(x = id_no,y = ESH_all,color = warmt)) +
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


# grid.draw(legend)


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
cowplot::ggsave(paste0(dir_mod,'figs/Figure_3_alt_meanensemble3.jpg'),fig,width = 220,height = 230,units='mm',dpi = 600)

#------------------------------------------------------------------------------
# check correlation body vs range size
d <- (data.frame(area=log10(tab$area),length=log10(tab$Length)))
cor(d[complete.cases(d),], method = "pearson")

summary(lm(log10(tab$area) ~ log10(tab$Length)))
library(ggplot2)
ggplot(tab,aes(x = log10(Length),y=log10(area))) +
  geom_point() +
  geom_smooth(method = 'lm')

summary(lm(log10(tab$area) ~ log10(tab$Length) + tab$code))
ggplot(tab,aes(x = log10(Length),y=log10(area),color = family)) +
  geom_point() +
  geom_smooth(method = 'lm')

library(lme4)
t<-lmer(log10(area)~log10(Length) + (1| order_/family), data= tab)
summary(t)

#-----------------------------------------------------------------------------
# check IUCN threat status boxplots

# this is available from 'Figure_SI_boxplot_scenario.R'
tab <- readRDS(paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds'))
tab <- droplevels( tab[tab$ESH_type == 'Total',] )

ncores = 8
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

# create table
tab <- merge(tab_median,iucn_simplified,by='id_no')

# IUCN code---
tab$code[tab$code == 'LR/lc'] <- 'LC'
tab$code[tab$code == 'EX'] <- 'DD'

tab$code <- factor(tab$code,levels=c('CR','EN','VU','LC','NT','DD'))
library(ggplot2)
ggplot(tab,aes(x = code, y = ESH_median,fill = warmt)) +
  geom_boxplot() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                     labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  ylab('% Range lost') +
  guides(fill = guide_legend(' ')) +
  xlab(' ') +
  theme_bw()
ggsave(paste0(dir_mod,'figs/Figure_SI_boxplot_IUCN_code.jpg'),width = 160,height = 160,units = 'mm',dpi=600)
