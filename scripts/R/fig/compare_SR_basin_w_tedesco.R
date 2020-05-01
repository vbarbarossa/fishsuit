library(sf); library(dplyr); library(ggplot2)

# get total number of species 
h <- readRDS('proc/ssp/points_template.rds')

# sample Tedesco basins
bas <- read_sf('data/Tedesco/Basin042017_3119.shp') %>% mutate(id = 1:nrow(.)) 
bas_ras <-  fasterize::fasterize(sf = bas,raster = raster(res = 1/12), field = 'id')
h$ws <- extract(bas_ras,h)
h <- h %>% filter(!is.na(ws))

# summarize species at Tedesco basin level the no species

# load species poly
sp <- read_sf('proc/species_ranges_merged.gpkg')

l <- st_intersects(h,sp,sparse = T)


t <- h %>%
  as_tibble() %>%
  select(-geometry) %>%
  mutate(row_no = 1:nrow(.))

tab <- lapply(unique(t$ws),
       function(x){
         tt <- t[t$ws == x,]
         return(
           data.frame(
             ws = x,
             no_sp = l[tt$row_no] %>% unlist %>% unique %>% length
           )
         )
       }) %>%
  do.call('rbind',.)

# merge with tedesco table
# use occurrences table
ted_occ <- read.csv('data/Tedesco/Occurrence_Table.csv',sep = ';') %>%
  as_tibble() %>%
  group_by('X1.Basin.Name') %>%
  summarise(no_sp_ted = length(unique(X5.Fishbase.Species.Code))) %>% 
  left_join(bas %>% select())





# hydrobasins level 8 within the Tedesco defined basins
hb8_ted <- readRDS('data/hb8_tedesco.rds')

# hydrobasins level 8 within each IUCN species
hb8_iucn <- lapply(list.files('ssp/ssp_hb8',full.names = T),function(x) readRDS(x)) %>%
  do.call('rbind',.)

# merge the two datasets, no need to keep species without a Tedesco basin
tab <- merge(hb8_iucn,hb8_ted,by='HYBAS_ID') # if too much merge them in the loop

# tabulate unique species within each basin
tab_unique <- foreach(bas = unique(tab$BasinName),.combine='rbind') %do% {
  t <- tab[tab$BasinName == bas,]
  
  cbind(t[1,which(colnames(t) == 'BasinName'):ncol(t)],
        data.frame(
          no_iucn = length(unique(t$binomial))
        )
  )
  
}

# merge tab_unique with species data on basins from Tedesco
ted_occ <- read.csv('data/Tedesco/Occurrence_Table.csv',sep=';')
# use fishbase valid species name -. 14953 species
ted_sp <- foreach(bas = unique(ted_occ$X1.Basin.Name),.combine='rbind') %do% data.frame(BasinName = bas, no_ted = length(unique(ted_occ[ted_occ$X1.Basin.Name == bas,'X5.Fishbase.Species.Code'])))

# merge with iucn occurrences to compare
mtab <- merge(tab_unique,ted_sp,by='BasinName')

# calculate difference

# load and merge with the shapefile for plotting
ws <- merge(read_sf('data/Tedesco/Basin042017_3119.shp'),mtab)
ws$ratio.perc <- ws$no_iucn/ws$no_ted*100
ws$diff.perc <- (ws$no_iucn - ws$no_ted)/ws$no_ted*100
ws$diff <- (ws$no_iucn - ws$no_ted)
plot(ws['ratio.perc'],breaks = c(seq(0,100,10)))
plot(ws[c('no_iucn','no_ted')])

ws$ratio.perc.corr <- ws$ratio.perc
ws$ratio.perc.corr[ws$ratio.perc > 100] <- 110
ws$ratio.perc.corr <- as.integer(round(ws$ratio.perc.corr,0))

ws.melt <- reshape2::melt(ws,measure.vars = c('no_iucn','no_ted'))
levels(ws.melt$variable) <- c('IUCN','Tedesco et al., 2017')


# base layers
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)


library(ggplot2); library(viridis)
p1 <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_sf(data = ws.melt%>% st_transform(crs_custom),aes(fill = value),alpha = 1,lwd = NA) +
  scale_fill_viridis_c(
    trans = 'log10',direction = -1
  ) +
  ggtitle('a') +
  facet_grid(variable~.) +
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

p2 <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  geom_sf(data = ws %>% st_transform(crs_custom),aes(fill = ratio.perc.corr),alpha = 1,lwd = NA) +
  scale_fill_viridis_c(
    breaks = seq(0,120,20),
    labels = c(seq(0,80,20),'>100',''),
    limits = c(0,120),
    discrete = F,option = 'C'
  ) +
  ggtitle('b') +
  facet_grid(variable~.) +
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

ggsave(paste0(dir_mod,'figs/maps_SR_compare_tedesco.jpg'),p1,
       width = 200,height = 200,dpi = 600,units = 'mm')
ggsave(paste0(dir_mod,'figs/maps_SR_compare_tedesco_diff.jpg'),p2,
       width = 200,height = 100,dpi = 600,units = 'mm')



p2 <- ggplot() +
  geom_sf(data = ws,aes(fill = ratio.perc.corr),alpha = 1) +
  scale_fill_viridis(breaks = seq(0,120,20),
                     labels = c(seq(0,80,20),'>100',''),
                     limits = c(0,120),expand = F,
                     discrete = F,option = 'C') +
  coord_sf(expand = F) +
  ggtitle('b') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(colour = 'transparent'),
        legend.position = 'bottom',
        # strip.background = element_rect('white'),
        legend.key.width = unit(8,'line'),
        legend.title = element_blank())

library(gridExtra)
p <- arrangeGrob(p1,p2,ncol = 1,nrow = 2,heights = c(2,1.2))
ggsave('D:/fishsuit/fishsuit_completeRun_warming_4targets/figs/Figure_SI_compare_SR_Tedesco.jpg',
       p,width = 220,height = 300,units = 'mm',dpi = 600)

saveRDS(ws,'data/compare_SR_tedesco.rds')
write.csv(as.data.frame(ws)[,-which(colnames(ws) == 'geometry')],'tabs/compare_SR_tedesco.csv',row.names = F)
