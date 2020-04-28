source('R/4targets/MASTER.R');

library(raster); library(foreach);


#> BASE DATA -------------------------------------------------------------------------------------
h <- readRDS('data/tab_template.rds')

# watersheds layer
ws_lyr <- raster('E:/GLOBIO-Aqua/output_5min_arcgis/basins_5min.tif')

# flow acc layer
fa_lyr <- raster('E:/GLOBIO-Aqua/output_5min_arcgis/flowAcc_5min.tif')

#area layer
area_lyr <- area(raster(res=1/12))

h@data$ws <- extract(ws_lyr,h)
h@data$fa <- extract(fa_lyr,h)
h@data$area <- extract(area_lyr,h)
#> ENSEMBLE SPECIES LOST -------------------------------------------------------------------------------------

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

dir_root <- 'fishsuit_completeRun_warming_4targets'

list.tabs.name <- paste0(dir_mod,dir_('figs/data/'),'tabs4zonal_stats.rds')

if(file.exists(list.tabs.name)){
  list.tabs <- readRDS(list.tabs.name)
}else{
  
  list.tabs <- foreach(n = seq_along(warming_targets)) %do% {
    
    atab <- attr.tab[attr.tab$warmt == warming_targets[n],]
    
    t <- as.data.frame(foreach(x = 1:nrow(atab),.combine='cbind') %do% {
      
      tab <- readRDS(paste0(dir_root,'/',atab$clmod[x],'/SR_tab_',
                            atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
      
      if(x == 1){
        return(
          cbind(tab@data$occ,(tab@data$occ - tab@data$all)/tab@data$occ)
        )
      }else{
        return(
          (tab@data$occ - tab@data$all)/tab@data$occ
        )
      }
    })
    colnames(t) <- c('occ',apply(atab,1,function(x) paste(x,collapse = '_')))
    tt <- t[-1]
    t[,paste0('me_',warming_targets[n])] <- apply(tt,1,function(x) median(x,na.rm=T))
    t[,paste0('av_',warming_targets[n])] <- apply(tt,1,function(x) mean(x,na.rm=T))
    t[,paste0('sd_',warming_targets[n])] <- apply(tt,1,function(x) sd(x,na.rm=T))
    
    return(t)
  }
  saveRDS(list.tabs,file=list.tabs.name)
}


#> ZONAL STATS ON WATERHSEDS -------------------------------------------------------------------------------------

# table of relative loss at different warming levels
tab <- cbind(data.frame(occ=list.tabs[[1]][,1]),
             do.call('cbind',lapply(list.tabs,function(x) x[,(ncol(x)-2):ncol(x)])),
             h@data[,-1])

by.ws <- split(h@data,f=h@data$ws)
# look up langerst watersheds
ws_size <- do.call('rbind',
                   lapply(by.ws,function(x) data.frame(ws_id = unique(x$ws),area = sum(x$area,na.rm=T)))
)

# and order based on size
ws_size <- ws_size[order(ws_size$area,decreasing = T),]

# # need to give a name to the watersheds
# rr <- ws_lyr
# rr[!rr %in% ws_size$ws_id[1:30]] <- NA
# plot(rr)

# # sample based on tedesco basins
library(sf)
wsted <- readRDS('data/compare_SR_tedesco.rds')
library(countrycode)
wsted$continent <- countrycode(sourcevar = wsted$Country,origin = 'country.name',destination = 'continent')
# # sample pcrglob basin ids 
# centr <- st_centroid(wsted)
# centr$ws_pcr <- extract(ws_lyr,centr)

# sample based on pcrglob basins
# ws_lyr_fil <- ws_lyr
# ws_lyr_fil[!ws_lyr_fil %in% ws_size$ws_id[1:500]] <- NA
# writeRaster(ws_lyr_fil,'figs/data/ws_lyr_filtered')
# poly <- rasterToPolygons(ws_lyr_fil,dissolve = T)
# rgdal::writeOGR(poly,'figs/data','ws_pcrglobwb_poly','ESRI Shapefile')

# centr2 <- st_centroid(read_sf('figs/data/ws_pcrglobwb_poly.shp'))
# library(fasterize)
# wsted$idted <- 1:nrow(wsted)
# c1 <- st_centroid(wsted)
# r <- fasterize(wsted,raster(res = 1/12),field = 'idted')
# 
# centr2$ws_ted <- extract(r,centr2)
# centr2$ws_ted <- as.integer(centr2$ws_ted)
# centr2_ <- merge(as.data.frame(centr2[!is.na(centr2$ws_ted),]),as.data.frame(wsted), by.x = 'ws_ted', by.y = 'idted')

centr2 <- st_centroid(wsted)

centr2$ws_id <- as.integer(raster::extract(ws_lyr,centr2))

# ws_size$names <- c(
#   'Amazon','Mississippi','Nile','Congo','Yenisey','Ob','Parana',
#   'Lena','Amur','Volga','Niger','Mackenzie','Yangtse','Ganges',
#   'Zambezi','Indus','Saint.Laurent','Lake.Eyre.Basin','Amu.Darya',
#   'Tchad','Orinoco','Tarim','Orange','Murray.Darling',
#   'Shatt.al.Arab','Yukon','Danube','Huanghe','Mekong','Tocantis',
#   (rep(NA,nrow(ws_size)-30))
# )

ws_size <- merge(ws_size,centr2,by = 'ws_id')
# filter out multiple records by selecting the largest ted watershed
ws_size <- foreach(i = unique(ws_size$ws_id),.combine = 'rbind') %do% {t = ws_size[ws_size$ws_id == i,]; return(t[which(t$Surf_area == max(t$Surf_area)),])}
# order to select largest watersheds
ws_size <- ws_size[order(ws_size$area,decreasing = T),]
row.names(ws_size) <- NULL

number_of_largest_basins = 200 #to consider
tab_ws <- foreach(i = ws_size$ws_id[1:number_of_largest_basins],.combine = 'rbind') %do% {
  
  t <- tab[tab$ws == i,]
  t <- t[!is.na(t$occ),]
  
  data.frame(ws_id = i,
             name = paste0(
               paste(do.call('c',strsplit(ws_size$BasinName[ws_size$ws_id == i],'\\.')),collapse = ' '),
               ' (',ws_size$no_iucn[ws_size$ws_id == i],'/',ws_size$no_ted[ws_size$ws_id == i],')'
               ),
             cover_ratio = ws_size$no_iucn[ws_size$ws_id == i]/ws_size$no_ted[ws_size$ws_id == i],
             richness_orig = ws_size$no_ted[ws_size$ws_id == i],
             area = ws_size$area[ws_size$ws_id == i],
             area_ted = ws_size$Surf_area[ws_size$ws_id == i],
             country = ws_size$Country[ws_size$ws_id == i],
             continent = ws_size$continent[ws_size$ws_id == i],
             ecoregion = ws_size$Ecoregion[ws_size$ws_id == i],
             endorheic = ws_size$Endorheic[ws_size$ws_id == i],
             av_1.5 = mean(t$me_1.5,na.rm=T),
             av_2.0 = mean(t$me_2.0,na.rm=T),
             av_3.2 = mean(t$me_3.2,na.rm=T),
             av_4.5 = mean(t$me_4.5,na.rm=T)
  )
}
str(tab_ws) #200 largest basins!

# order based on original richness according to tedesco et al.
tab_ws <- tab_ws[order(tab_ws$richness_orig,decreasing = T),]

tab_ws_total <- tab_ws

# filter out watersheds with a coverage lower than 0.5
# and select the 5 richest per continent
tab_ws <- tab_ws[tab_ws$cover_ratio>=0.5 | tab_ws$name == 'Amazon (404/2273)',] #include Amazon anyway
tab_ws <- do.call('rbind',lapply(split(tab_ws,tab_ws$continent),function(x) x[x$richness_orig >= 50,][1:6,]))
tab_ws <- tab_ws[!is.na(tab_ws$ws_id),]

# save pre-processed watershed data for alternative approach
saveRDS(by.ws[sort(tab_ws_total$ws_id)],'data/preprocess_altFig2/splitted_watersheds200.rds')
write.csv(tab_ws_total,'data/preprocess_altFig2/tab_watersheds200.csv',row.names = F)

# barplots --------------------------------------------------------------------------------------------------------
# only relative

# for main figure (tab_ws)

library(ggplot2); library(reshape2); library(RColorBrewer)
df <- melt(tab_ws[,-which(colnames(tab_ws) %in% c('cover_ratio','richness_orig'))],id.vars = c('ws_id','name','area','area_ted','country','ecoregion','endorheic','continent'))
levels(df$variable) <- c('1.5C','2.0C','3.2C','4.5C')
df <- df[order(df$value, decreasing = T),]

# order names from highest to lowest impact
df$name <- factor(df$name,levels = rev(as.character(df[df$variable == '4.5C','name'])))

# levels(df$ecoregion) <- c(levels(df$ecoregion)[-6],'Au.')
# levels(df$continent) <- c(levels(df$continent)[-5],'Oc.')

p <- ggplot(df,aes(y = value, x = name, fill = variable)) +
  geom_bar(stat='identity',position = 'identity',width = 0.5) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                    labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  scale_y_continuous(breaks = seq(0,1,0.25),labels = seq(0,1,0.25)) +
  guides(fill = guide_legend(title=NULL)) +
  ylab('Average Local Cumulative Range Loss [-]') +
  xlab(' ') +
  coord_flip(ylim = c(0,1.0001),expand = F) +
  facet_grid(continent~.,scales = 'free_y',space = 'free_y',switch = 'y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        legend.position = c(0.92,0.12),
        strip.placement = 'outside',
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"))
p
ggsave(paste0(dir_mod,'figs/Figure_2_richest_basins_per_continent_alt.jpg'),p,width = 160,height = 140,units='mm',dpi = 600)


# and the total ----
df <- melt(tab_ws_total[,-which(colnames(tab_ws) %in% c('cover_ratio','richness_orig'))],id.vars = c('ws_id','name','area','area_ted','country','ecoregion','endorheic','continent'))
levels(df$variable) <- c('1.5C','2.0C','3.2C','4.5C')
df <- df[order(df$value, decreasing = T),]

# order names from highest to lowest impact
df$name <- factor(df$name,levels = rev(as.character(df[df$variable == '4.5C','name'])))

p <- ggplot(df,aes(y = value, x = name, fill = variable)) +
  geom_bar(stat='identity',position = 'identity',width = 0.5) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                    labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  scale_y_continuous(breaks = seq(0,1,0.25),labels = c('',seq(0.25,1,0.25))) +
  guides(fill = guide_legend(title=NULL)) +
  ylab('Average Local Cumulative Range Loss [-]') +
  xlab(' ') +
  coord_flip(ylim = c(0,1.0001),expand = F) +
  facet_grid(continent~.,scales = 'free_y',space = 'free_y',switch = 'y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        legend.position = c(0.92,0.03),
        strip.placement = 'outside',
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"))
# p
ggsave(paste0(dir_mod,'figs/Figure_SI_richest_basins_per_continent_200alt2.jpg'),p,width = 200,height = 600,units='mm',dpi = 600)
ggsave(paste0(dir_mod,'figs/Figure_SI_richest_basins_per_continent_200t2.jpg'),p,width = 100,height = 600,units='mm',dpi = 600)


# absolute also NOT IMPLEMENTED FOR NOW

tab_ws <- foreach(i = ws_size$ws_id[1:number_of_largest_basins],.combine = 'rbind') %do% {
  
  t <- tab[tab$ws == i,]
  t <- t[!is.na(t$occ),]
  rbind(
    data.frame(ws_id = i,
               name = paste0(
                 paste(do.call('c',strsplit(ws_size$BasinName[ws_size$ws_id == i],'\\.')),collapse = ' '),
                 ' (',ws_size$no_iucn[ws_size$ws_id == i],'/',ws_size$no_ted[ws_size$ws_id == i],')'
               ),
               area = ws_size$area[ws_size$ws_id == i],
               area_ted = ws_size$Surf_area[ws_size$ws_id == i],
               country = ws_size$Country[ws_size$ws_id == i],
               continent = ws_size$continent[ws_size$ws_id == i],
               ecoregion = ws_size$Ecoregion[ws_size$ws_id == i],
               endorheic = ws_size$Endorheic[ws_size$ws_id == i],
               av_1.5 = mean(t$me_1.5,na.rm=T),
               av_2.0 = mean(t$me_2.0,na.rm=T),
               av_3.2 = mean(t$me_3.2,na.rm=T),
               av_4.5 = mean(t$me_4.5,na.rm=T),
               type = 'Relative'
    ),
    data.frame(ws_id = i,
               name = paste0(
                 paste(do.call('c',strsplit(ws_size$BasinName[ws_size$ws_id == i],'\\.')),collapse = ' '),
                 ' (',ws_size$no_iucn[ws_size$ws_id == i],'/',ws_size$no_ted[ws_size$ws_id == i],')'
               ),
               area = ws_size$area[ws_size$ws_id == i],
               area_ted = ws_size$Surf_area[ws_size$ws_id == i],
               country = ws_size$Country[ws_size$ws_id == i],
               continent = ws_size$continent[ws_size$ws_id == i],
               ecoregion = ws_size$Ecoregion[ws_size$ws_id == i],
               endorheic = ws_size$Endorheic[ws_size$ws_id == i],
               av_1.5 = mean(t$me_1.5*t$occ,na.rm=T),
               av_2.0 = mean(t$me_2.0*t$occ,na.rm=T),
               av_3.2 = mean(t$me_3.2*t$occ,na.rm=T),
               av_4.5 = mean(t$me_4.5*t$occ,na.rm=T),
               type = 'Absolute'
    )
  )
  
}

library(ggplot2); library(reshape2); library(RColorBrewer)
df <- melt(tab_ws[-c(3,4,5,7,8)],id.vars = c('ws_id','name','type','continent'))
levels(df$variable) <- c('1.5C','2.0C','3.2C','4.5C')
df <- df[order(df$value, decreasing = T),]

# levels(df$ecoregion) <- c(levels(df$ecoregion)[-6],'Au.')
levels(df$continent) <- c(levels(df$continent)[-5],'Oc.')


# order names from highest to lowest impact
df$name <- factor(df$name,levels = rev(as.character(df[df$variable == '4.5C' & df$type == 'Relative','name'])))

p <- ggplot(df,aes(y = value, x = name, fill = variable)) +
  geom_bar(stat='identity',position = 'identity',width = 0.5) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                    labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  # scale_y_continuous(breaks = seq(0,1,0.25),labels = c('',seq(0.25,1,0.25))) +
  guides(fill = guide_legend(title=NULL)) +
  xlab(' ') +
  ylab(' ') +
  coord_flip(expand = F) +
  facet_grid(continent~type,scales = 'free',space = 'free_y',switch = 'both') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        legend.position = 'none', #c(0.85,0.5)
        strip.placement = 'outside',
        strip.background = element_blank(),
        panel.spacing = unit(1.3, "lines"))


ggsave(paste0(dir_mod,'figs/Figure_3_major_basins2_cont_abs.jpg'),p,width = 180,height = 250,units='mm',dpi = 600)








#> ZONAL STATS ON FLOWACC -------------------------------------------------------------------------------------

#need to find a way to cluster the upstream area
library(ggplot2); library(reshape2); library(RColorBrewer)
df <- melt(tab[!is.na(tab$occ),c('occ','fa',paste0('me_',c('1.5','2.0','3.2','4.5')))],id.vars = c('occ','fa'))
df$fa <- df$fa+1

df$fa_cat <- 0
for(i in seq(0.5,4.5,1)){
  df$fa_cat[log10(df$fa) >= i & log10(df$fa) < (i+0.5)] <- (i+0.5)
}
df$fa_cat <- factor(df$fa_cat)

#violin is ugly
#boxplot hmmm
ggplot(df,aes(x = fa_cat,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                    labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  guides(fill = guide_legend(title=NULL)) +
  xlab('LOG-upstream grid cells') +
  ylab('Relative species lost') +
  theme_bw()

# no trend
summary(lm(value ~ log10(fa), data = df))

#scatter
p <- ggplot(df[df$variable == 'me_3.2',],aes(x = log10(fa),y = value)) +
  geom_point(alpha = 0.1,size = 0.1,color = 'Grey') +
  xlab('LOG-upstream number of grid cells') +
  ylab('Relative species lost') +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(paste0(dir_mod,'figs/Figure_SI_flowacc_vs_relative_lost.jpg'),p,width = 120,height = 120,units='mm',dpi = 600)


#> ZONAL STATS ON AV.FLOW -------------------------------------------------------------------------------------
# based on qav historical

#average historical flow raster layers
qav_med <- calc(brick(foreach(cl = climate_models) %do% raster(paste0(dir_mod,cl,'/pcrglobwb_processed/merged/Qav_hist.tif'))),
                fun = median,na.rm = T)

# extract the values
h@data$qav <- extract(qav_med,h)

tab <- cbind(data.frame(occ=list.tabs[[1]][,1]),
             do.call('cbind',lapply(list.tabs,function(x) x[,(ncol(x)-2):ncol(x)])),
             h@data[,-1])

df <- melt(tab[!is.na(tab$occ),c('occ','qav',paste0('me_',c('1.5','2.0','3.2','4.5')))],id.vars = c('occ','qav'))
df <- df[df$qav > 0,] # remove zero-flow cells
df$log_qav <- log10(df$qav)

# no correlation
summary(lm(value ~ log_qav, data = df[df$variable == 'me_1.5',]))


p <- ggplot(df,aes(x = log_qav,y = value)) +
  geom_point(alpha = 0.1,size = 0.1,color = 'Grey') +
  xlab('LOG-average annual flow (median over historical runs)') +
  ylab('Relative species lost') +
  theme_bw() +
  facet_wrap('variable',ncol=2,nrow=2) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        strip.background = element_blank())

ggsave(paste0(dir_mod,'figs/Figure_SI_qav_vs_relative_lost.jpg'),p,width = 220,height = 220,units='mm',dpi = 600)





