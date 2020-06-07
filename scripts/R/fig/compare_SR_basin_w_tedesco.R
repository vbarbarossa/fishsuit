library(sf); library(dplyr); library(ggplot2); library(raster)

hbas <- read_sf('data/hybas12_points_nolakes.gpkg')

# sample Tedesco basins
bas <- read_sf('data/Tedesco/Basin042017_3119.shp') %>% mutate(id = 1:nrow(.)) 
bas_ras <-  fasterize::fasterize(sf = bas,raster = raster(res = 1/12), field = 'id')
hbas$ws <- extract(bas_ras,hbas)
hbas_tab <- hbas %>% as_tibble() %>% dplyr::select(-geom) %>% filter(!is.na(ws))

# summarize species at Tedesco basin level the no species
sp <- readRDS('proc/species_ranges_raw_on_hybas12.rds')

tab <- left_join(sp,hbas_tab %>% dplyr::select(HYBAS_ID,ws)) %>%
  group_by(ws) %>%
  summarize(sp_no = length(unique(binomial)))

tab <- left_join(tab,bas %>% as_tibble() %>% dplyr::select(BasinName,ws = id,-geometry))

# merge with tedesco table
# use occurrences table
ted_occ <- read.csv('data/Tedesco/Occurrence_Table.csv',sep = ';') %>%
  as_tibble() %>%
  rename(BasinName = X1.Basin.Name) %>%
  group_by(BasinName) %>%
  summarise(no_sp_ted = length(unique(X5.Fishbase.Species.Code))) %>% 
  left_join(bas %>% as_tibble() %>% dplyr::select(BasinName,ws = id,-geometry))

tab <- inner_join(tab,ted_occ)

# calculate difference

# load and merge with the shapefile for plotting
ws <- inner_join(bas,tab) %>%
  mutate(
    ratio.perc = sp_no/no_sp_ted * 100,
    diff.perc = (sp_no - no_sp_ted)/no_sp_ted * 100,
    diff = sp_no - no_sp_ted
  )

ws$ratio.perc.corr <- ws$ratio.perc
ws$ratio.perc.corr[ws$ratio.perc > 100] <- 110
ws$ratio.perc.corr <- as.integer(round(ws$ratio.perc.corr,0))

ws.melt <- reshape2::melt(ws %>% as_tibble() %>% dplyr::select(-geometry),measure.vars = c('sp_no','no_sp_ted'))
levels(ws.melt$variable) <- c('This study','Tedesco et al., 2017')

ws.melt <- left_join(bas,ws.melt,by = 'BasinName') %>% filter(!is.na(value))
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
  geom_sf(data = ws.melt %>% st_transform(crs_custom),aes(fill = value),alpha = 1,lwd = NA) +
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
    breaks = c(seq(0,100,20),110),
    labels = c(seq(0,80,20),'>100',''),
    limits = c(0,110),
    option = 'C'
  ) +
  ggtitle('b') +
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

ggsave(paste0('figs/maps_SR_compare_tedesco.jpg'),p1,
       width = 200,height = 200,dpi = 600,units = 'mm')
ggsave(paste0('figs/maps_SR_compare_tedesco_diff.jpg'),p2,
       width = 200,height = 100,dpi = 600,units = 'mm')


library(gridExtra)
p <- arrangeGrob(p1,p2,ncol = 1,nrow = 2,heights = c(2,1.2))
ggsave('figs/compare_SR_Tedesco_combined.jpg',
       p,width = 220,height = 300,units = 'mm',dpi = 600)

saveRDS(ws,'proc/compare_SR_tedesco.rds')
write.csv(as.data.frame(ws)[,-which(colnames(ws) == 'geometry')],'proc/compare_SR_tedesco.csv',row.names = F)
