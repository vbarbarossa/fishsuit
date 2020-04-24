source('config_local.R'); # always load as functions are loaded within this script

libinv(c('dplyr','sf'))

# species ranges
ranges <- read_sf('proc/species_ranges_merged.gpkg') %>%
  as_tibble() %>% select(-geom)

# # # # load species habitat type <<<<<<< need to check why we are losing almost 100 species when merging
hab <- left_join(ranges,read.csv('proc/species_traits.csv'))

# override vars
vars <- c('Qmi','Qma','Qzf',
          'Qcv','Qcv',
          # 'Qve',
          'Tma','Tmi',
          'Tcv','Tcv')
thresholds <- c('2.5','100.0','97.5',
                '2.5','97.5',
                # '100.0',
                '97.5','2.5',
                '2.5','97.5')

# make one table per GCM and then average across them (include only common ids)
niche_list <- lapply(
  seq_along(vars),
  function(j){
    
    t <- lapply(
      seq_along(climate_models),
      function(x){
        
        d <- readRDS(paste0('proc/',climate_models[x],'/niches/',vars[j],'.rds')) %>% 
          select(ID,paste0(thresholds[j],'%')) %>% as_tibble()
        if(x == 1){return(d)}else{return(d[,2])}
      }) %>% do.call('cbind',.)
    
    res <- data.frame(id_no = t$ID)
    res[,paste0(vars[j],'_',thresholds[j],'%')] <- apply(t[-1],1,function(x) mean(x,na.rm=T))
    
    return(res)
  }
)

# make sure all ids correspond properly
niche <- niche_list[[1]] %>% as_tibble()
for(i in seq_along(vars)[-1]){
  niche <- inner_join(niche,niche_list[[i]])
}

# join with habitat table
niche <- left_join(niche,hab %>% mutate(id_no = as.character(id_no)))

# join with number of grid cells
niche <- left_join(niche,readRDS(paste0('proc/',climate_models[1],'/niches/','Qmi.rds')) %>% select(id_no = ID, no_cells = no.grids))

# save table
write.csv(niche,'proc/thresholds_average_all.csv',row.names = F)

#---------------------------------------------------------------------------------------------
# FILTERING

nrow(niche)
# [1] 12934

# 1 - filter out exclusively lentic species

niche <- niche %>%
  filter(lentic_only == 0)
  
nrow(niche)  
# [1] 11774

# 2 - min no grid cells
niche <- niche %>%
  filter(no_cells >= min_no_grid_cells,
         !is.na(`Qmi_2.5%`),!is.nan(`Qmi_2.5%`),!is.infinite(`Qmi_2.5%`))

nrow(niche)
# [1] 11425

write.csv(niche,'proc/thresholds_average_filtered.csv',row.names = F)

#---------------------------------------------------------------------------------------------
# check correlations of thresholds

libinv(c('ggplot2','purrr','tidyr','bestNormalize','foreach'))

tab <- niche[,1:(length(vars)+1)]

apply(tab %>% select(-id_no),2,function(x) sum(is.na(x)))
apply(tab %>% select(-id_no),2,function(x) sum(x == 0))

# seems that some Tmi value are set to 0, then set them to the minimum 273.15
tab$`Tmi_2.5%`[tab$`Tmi_2.5%` < 275.15] <- 273.15

# check distribution
(hist <- tab[-1] %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal())


# need to transform variables to correctly estimate correlation metrics
BN <- list()
tabN <- foreach(i = seq_along(vars),.combine = 'cbind') %do% {
  
  BN[[i]] <- list()
  BN[[i]][[1]] <- bestNormalize(as.data.frame(tab)[,i+1],allow_orderNorm = FALSE)
  
  val <- scale(BN[[i]][[1]]$x.t)
  BN[[i]][[2]] <- attr(val,"scaled:center")
  BN[[i]][[3]] <- attr(val,"scaled:scale")
  
  names(BN[[i]]) <- c('BN','center','scale')
  
  d <- as.data.frame(as.numeric(val))
  colnames(d) <- colnames(tab)[i+1] 
  
  return( d )
  
}
names(BN) <- colnames(tab)[-1]

# and check again the distribution
(histN <- tabN %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal())


#' First bivariate correlations
cmN <- cor(tabN, use = "pairwise.complete.obs")
# and visualize
corrplot::corrplot(cmN, method = 'number', type = 'lower', number.cex = 1)

#' First bivariate correlations
cm <- cor(tab[-1], use = "pairwise.complete.obs")
# and visualize
corrplot::corrplot(cm, method = 'number', type = 'lower', number.cex = 1)

# make sure a figs directory exists
dir_('figs')

# save figures

# histograms
ggsave('figs/hist_covariates.jpg',hist,width = 160, height = 160,dpi = 600, units = 'mm')

ggsave('figs/hist_covariates_normalized.jpg',histN,width = 160, height = 160,dpi = 600, units = 'mm')

# corrplots
jpeg('figs/corrplot_covariates.jpg',width = 160, height = 160, res = 600, units = 'mm')
corrplot::corrplot(cm, method = 'number', type = 'lower', number.cex = 1)
dev.off()

jpeg('figs/corrplot_covariates_normalized.jpg',width = 160, height = 160, res = 600, units = 'mm')
corrplot::corrplot(cmN, method = 'number', type = 'lower', number.cex = 1)
dev.off()
