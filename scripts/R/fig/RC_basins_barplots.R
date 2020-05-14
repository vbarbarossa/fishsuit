source('config.R'); # always load as functions are loaded within this script

library(raster); library(foreach); library(dplyr); library(sf)


#> BASE DATA -------------------------------------------------------------------------------------
h <- readRDS('proc/ssp/points_template.rds')

# watersheds layer
ws_lyr <- raster('data/pcrglobwb_hydrography/basins_5min.tif')

h$ws <- extract(ws_lyr,h)
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

lst_tabs <- foreach(n = seq_along(warming_targets)) %do% {
  
  atab <- attr.tab[attr.tab$warmt == warming_targets[n],]
  
  t <- as.data.frame(foreach(x = 1:nrow(atab),.combine='cbind') %do% {
    
    tab <- readRDS(paste0('proc/',atab$clmod[x],'/SR_tab_',
                          atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
    
    if(x == 1){
      return(
        cbind(tab$occ,(tab$occ - tab$all)/tab$occ)
      )
    }else{
      return(
        (tab$occ - tab$all)/tab$occ
      )
    }
  })
  colnames(t) <- c('occ',apply(atab,1,function(x) paste(x,collapse = '_')))
  
  t[,paste0('me_',warming_targets[n])] <- matrixStats::rowMedians(as.matrix(t[-1]),na.rm = T)
  t[,paste0('av_',warming_targets[n])] <- rowMeans(t[-1],na.rm=T)
  
  return(t)
}

lst_tabs_dsp <- foreach(n = seq_along(warming_targets)) %do% {
  
  atab <- attr.tab[attr.tab$warmt == warming_targets[n],]
  
  t <- as.data.frame(foreach(x = 1:nrow(atab),.combine='cbind') %do% {
    
    tab <- readRDS(paste0('proc/',atab$clmod[x],'/SR_tab_dispersal_',
                          atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
    
    if(x == 1){
      return(
        cbind(tab$occ,(tab$occ - tab$all)/tab$occ)
      )
    }else{
      return(
        (tab$occ - tab$all)/tab$occ
      )
    }
  })
  colnames(t) <- c('occ',apply(atab,1,function(x) paste(x,collapse = '_')))
  
  t[,paste0('me_',warming_targets[n])] <- matrixStats::rowMedians(as.matrix(t[-1]),na.rm = T)
  t[,paste0('av_',warming_targets[n])] <- rowMeans(t[-1],na.rm=T)
  
  return(t)
}

#> ZONAL STATS ON WATERHSEDS -------------------------------------------------------------------------------------

# table of relative loss at different warming levels
tab <- cbind(data.frame(occ=lst_tabs[[1]][,1]),
             do.call('cbind',lapply(lst_tabs,function(x) x[,(ncol(x)-1):ncol(x)])),
             h %>% as_tibble() %>% dplyr::select(-geometry))
tab_dsp <- cbind(data.frame(occ=lst_tabs_dsp[[1]][,1]),
                 do.call('cbind',lapply(lst_tabs_dsp,function(x) x[,(ncol(x)-1):ncol(x)])),
                 h %>% as_tibble() %>% dplyr::select(-geometry))


by.ws <- split(h,f=h$ws)
# look up langerst watersheds
ws_size <- do.call(
  'rbind',
  lapply(by.ws,function(x) data.frame(ws_id = unique(x$ws),area = sum(x$area,na.rm=T)))
)

# and order based on size
ws_size <- ws_size[order(ws_size$area,decreasing = T),]

# # sample based on tedesco basins
wsted <- readRDS('proc/compare_SR_tedesco.rds')
library(countrycode)
wsted$continent <- countrycode(sourcevar = wsted$Country,origin = 'country.name',destination = 'continent')

centr2 <- st_centroid(wsted)

centr2$ws_id <- as.integer(raster::extract(ws_lyr,centr2))

ws_size <- merge(ws_size,centr2,by = 'ws_id')
# filter out multiple records by selecting the largest ted watershed
ws_size <- foreach(i = unique(ws_size$ws_id),.combine = 'rbind') %do% {t = ws_size[ws_size$ws_id == i,]; return(t[which(t$Surf_area == max(t$Surf_area)),])}
# order to select largest watersheds
ws_size <- ws_size[order(ws_size$area,decreasing = T),]
row.names(ws_size) <- NULL

number_of_largest_basins = 200 #to consider
d1 <- foreach(i = ws_size$ws_id[1:number_of_largest_basins],.combine = 'rbind') %do% {
  
  t <- tab[tab$ws == i,]
  t <- t[!is.na(t$occ),]
  
  
  data.frame(ws_id = i,
             name = paste0(
               paste(do.call('c',strsplit(ws_size$BasinName[ws_size$ws_id == i],'\\.')),collapse = ' '),
               ' (',ws_size$sp_no[ws_size$ws_id == i],')'
             ),
             cover_ratio = ws_size$sp_no[ws_size$ws_id == i]/ws_size$no_sp_ted[ws_size$ws_id == i],
             richness_orig = ws_size$no_sp_ted[ws_size$ws_id == i],
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
             scenario = 'no dispersal'
  )
}

# order based on original richness according to tedesco et al.
d1 <- d1[order(d1$richness_orig,decreasing = T),]

d1_tot <- d1

# filter out watersheds with a coverage lower than 0.5
# and select the 5 richest per continent
d1 <- do.call('rbind',lapply(split(d1,d1$continent),function(x) x[x$richness_orig >= 50,][1:6,]))
d1 <- d1[!is.na(d1$ws_id),]

d2 <- foreach(i = ws_size$ws_id[1:number_of_largest_basins],.combine = 'rbind') %do% {
  
  
  t <- tab_dsp[tab_dsp$ws == i,]
  t <- t[!is.na(t$occ),]
  
  
  data.frame(ws_id = i,
             name = paste0(
               paste(do.call('c',strsplit(ws_size$BasinName[ws_size$ws_id == i],'\\.')),collapse = ' '),
               ' (',ws_size$sp_no[ws_size$ws_id == i],')'
             ),
             cover_ratio = ws_size$sp_no[ws_size$ws_id == i]/ws_size$no_sp_ted[ws_size$ws_id == i],
             richness_orig = ws_size$no_sp_ted[ws_size$ws_id == i],
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
             scenario = 'maximal dispersal'
  )
}

# order based on original richness according to tedesco et al.
d2 <- d2[order(d2$richness_orig,decreasing = T),]

d2_tot <- d2

# filter out watersheds with a coverage lower than 0.5
# and select the 5 richest per continent
d2 <- do.call('rbind',lapply(split(d2,d2$continent),function(x) x[x$richness_orig >= 50,][1:6,]))
d2 <- d2[!is.na(d2$ws_id),]

ws_tab_tot <- rbind(d1_tot,d2_tot)
ws_tab <- rbind(d1,d2)


# save pre-processed watershed data for alternative approach
write.csv(ws_tab_tot %>% dplyr::select(name,continent,av_1.5,av_2.0,av_3.2,av_4.5,scenario),paste0('figshare/RC_basins.csv'),row.names = F)

# barplots --------------------------------------------------------------------------------------------------------
# only relative

# for main figure (ws_tab)

library(ggplot2); library(reshape2); library(RColorBrewer)
df <- melt(ws_tab %>% dplyr::select(-c('cover_ratio','richness_orig')),measure.vars = c('av_1.5','av_2.0','av_3.2','av_4.5'))
levels(df$variable) <- c('1.5C','2.0C','3.2C','4.5C')
df <- df[order(df$value, decreasing = T),]
df <- droplevels(df)

# order names from highest to lowest impact
df$name <- factor(df$name,levels = rev(unique(as.character(df[df$variable == '4.5C','name']))))

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
  facet_grid(continent ~ scenario,scales = 'free_y',space = 'free_y',switch = 'y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        # legend.direction = 'horizontal',
        legend.position = c(0.92,0.12),
        strip.placement = 'outside',
        strip.background = element_blank(),
        panel.spacing = unit(1.2, "lines"))
p
pg <- ggplotGrob(p)

for(i in which(grepl("strip", pg$layout$name))){
  pg$grobs[[i]]$layout$clip <- "off"
}

ggsave(paste0('figs/barplot_basins_selection.jpg'),pg,width = 183,height = 150,units='mm',dpi = 600)
ggsave(paste0('figs/barplot_basins_selection.pdf'),pg,width = 183,height = 150,units='mm')

# and the total ----
df <- melt(ws_tab_tot %>% dplyr::select(-c('cover_ratio','richness_orig')),measure.vars = c('av_1.5','av_2.0','av_3.2','av_4.5'))
levels(df$variable) <- c('1.5C','2.0C','3.2C','4.5C')
df <- droplevels(df)
df <- df[order(df$value, decreasing = T),]

# order names from highest to lowest impact
df$name <- factor(df$name,levels = rev(unique(as.character(df[df$variable == '4.5C','name']))))

p <- ggplot(df,aes(y = value, x = name, fill = variable)) +
  geom_bar(stat='identity',position = 'identity',width = 0.5) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                    labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  scale_y_continuous(breaks = seq(0,1,0.25),labels = c('',seq(0.25,1,0.25))) +
  guides(fill = guide_legend(title=NULL)) +
  ylab('Average Local Cumulative Range Loss [-]') +
  xlab(' ') +
  coord_flip(ylim = c(0,1.0001),expand = F) +
  facet_grid(continent~scenario,scales = 'free_y',space = 'free_y',switch = 'y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        legend.position = c(0.92,0.027),
        strip.placement = 'outside',
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"))
# p
ggsave(paste0('figs/barplot_basins.jpg'),p,width = 200,height = 600,units='mm',dpi = 1000)
ggsave(paste0('figs/barplot_basins.pdf'),p,width = 200,height = 600,units='mm')



