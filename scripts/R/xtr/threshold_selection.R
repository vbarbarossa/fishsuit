source('config_local.R'); # always load as functions are loaded within this script

libinv(c('dplyr'))

# species ranges
ranges <- read_sf('proc/species_ranges_merged.gpkg') %>%
  as_tibble() %>% select(-geom)

# # # # load species habitat type <<<<<<< need to check why we are losing almost 100 species when merging
hab <- left_join(ranges,read.csv('proc/species_traits.csv'))

for(g in seq_along(climate_models)){
  
  # override vars
  vars <- c('Qmi','Qma','Qzf','Tma','Tmi')
  thresholds <- c('2.5','100.0','97.5','97.5','2.5')
  clmod = climate_models[g]
  
  dir_niche <- paste0('proc/',clmod,'/niches/')
  
  
  niche_list <- list()
  for(j in seq_along(vars)){
    niche <- readRDS(paste0(dir_niche,vars[j],'.rds'))
    
    niche <- niche[,c('ID',paste0(thresholds[j],'%'))]
    colnames(niche)[2] <- paste0(vars[j],'_',colnames(niche)[2])
    
    niche_list[[j]] <- niche
  }  
  niche <- as.data.frame(do.call('cbind',
                                 lapply(
                                   niche_list,function(x) x[,2]
                                 )
  ))
  colnames(niche) <- sapply(niche_list,function(x) colnames(x)[2])
  niche <- cbind(niche_list[[1]][,1],niche)
  colnames(niche)[1] <- 'id_no'
  
  #---------------------------------------------------------------------------------------------
  # FILTERING
  # save ids step by step
  cat(paste0(rep('-',30)),'\n')
  cat('FILTERING:\n')
  
  # STEP 0 ###
  cat('STEP 0: initial species no. = ',nrow(niche),'\n')
  
  # STEP 1 ###
  # filter out species outside the pcrglobwb domain and
  # with less than a user defined number of grid-cells threshold
  d <- readRDS(paste0(dir_niche,'Qmi.rds'))
  step1 <- d[!is.na(d$mean) & d$no.grids >= min_no_grid_cells,1]
  
  cat('STEP 1: removing ',(nrow(niche) - length(step1)),' species represented by less than ',min_no_grid_cells,' grid-cell or falling out of the pcr-globwb domain\n')
  niche <- niche[niche$id_no %in% step1,]
  
  
  # STEP 2
  ## merge habitat
  niche <- left_join(niche %>% as_tibble,hab %>% mutate(id_no = as.character(id_no)))
  
  # STEP 3
  ## filter for lotic sp only
  niches_f <- niche[niche$lentic_only == 0,]
  cat('removing ',(nrow(niche) - nrow(niches_f)),' species exclusively lentic\n')
  
  niches_f <- niches_f[niches_f$`Tmi_2.5%` > 273.15,] #<<<<<<<<<< WHAT? rather set it to this value where it is lower?
  
  write.csv(niches_f,paste0(dir_niche,'niches_filtered.csv'),row.names = F)
  
  cat(paste0(rep('-',30)),'\n\n')
  
  
}

libinv(c('ggplot2','purrr','tidyr','bestNormalize','foreach'))

tab <- niches_f %>%
  as_tibble() %>%
  .[1:5]

apply(tab %>% select(-id_no),2,function(x) sum(is.na(x)))
apply(tab %>% select(-id_no),2,function(x) sum(x == 0))


# check distribution

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
tabN %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_bw()


#' First bivariate correlations
cm <- cor(tabN, use = "pairwise.complete.obs")
# and visualize
corrplot::corrplot(cm, method = 'number', type = 'lower', number.cex = 1)

with(tab,cor(`Qzf_97.5%`,`Qmi_2.5%`))

