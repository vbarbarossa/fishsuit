source('config_local.R'); # always load as functions are loaded within this script

libinv(c('dplyr'))

for(g in seq_along(climate_models)){
  
  # override vars
  vars <- c('Qmi','Qma','Qzf','Tma','Tmi')
  thresholds <- c('2.5','100.0','97.5','97.5','2.5')
  clmod = climate_models[g]
  scen = scenarios[scenarios != 'hist'] #<
  years = floor(mean(timespan_scen)) #<
  
  dir_niche <- paste0(dir_model,'proc/',clmod,'/niches/')
  
  # # # # load iucn data
  iucn <- lapply(1:2, function(i) foreign::read.dbf(paste0('data/FW_FISH_20181113/FW_FISH_PART_',i,'.dbf'))) %>%
    do.call('rbind',.)
  iucn <- iucn[,c('id_no','binomial')][!duplicated(iucn[,c('id_no','binomial')]),]
  row.names(iucn) <- 1:nrow(iucn)
  
  # # # # load species habitat type <<<<<<< need to check why we are losing almost 100 species when merging
  hab <- merge(iucn,read.csv('iucn_habitat_type.csv'),by='binomial')
  
  niche_list <- list()
  for(j in seq_along(vars)){
    niche <- readRDS(paste0(dir_niche,vars[j],'.rds'))
      
    
    # remove ranges with NAs
    # remove ranges with SD = 0, as they are sampled for only one 5 arcmin cell
    # niche <- niche[!is.na(niche$mean),] # & niche$sd != 0
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
  colnames(niche)[1] <- 'IUCN_ID'
  
  #---------------------------------------------------------------------------------------------
  # FILTERING
  # save ids step by step
  cat(paste0(rep('-',30)),'\n')
  cat('FILTERING:\n')
  
  # STEP 0 ###
  cat('STEP 0: initial species no. = ',nrow(niche),'\n')
  saveRDS(niche$IUCN_ID,paste0(dir_niche,'filtering_step0_IDs.rds'))
  
  # STEP 1 ###
  # filter out species outside the pcrglobwb domain and
  # with less than a user defined number of grid-cells threshold
  d <- readRDS(paste0(dir_niche,'Qmi.rds'))
  step1 <- d[!is.na(d$mean) & d$no.grids >= min_no_grid_cells,1]
  
  cat('STEP 1: removing ',(nrow(niche) - length(step1)),' species represented by less than ',min_no_grid_cells,' grid-cell or falling out of the pcr-globwb domain\n')
  niche <- niche[niche$IUCN_ID %in% step1,]
  saveRDS(niche$IUCN_ID,paste0(dir_niche,'filtering_step1_IDs.rds'))
  
  niches_f <- niche
  
  if(filter_lentic_out){
    
    # STEP 2
    ## merge habitat
    niches <- merge(niche,hab[,c('id_no','binomial','lotic','lentic')],by.x='IUCN_ID',by.y='id_no')
    cat('removing ',(nrow(niche) - nrow(niches)),' species when merging with habitat info\n')
    saveRDS(niches$IUCN_ID,paste0(dir_niche,'/filtering_step2_IDs.rds'))
    
    # STEP 3
    ## filter for lotic sp only
    niches_f <- niches[(niches$lotic == 1),]
    cat('removing ',(nrow(niches) - nrow(niches_f)),' species exclusively lentic\n')
    saveRDS(niches_f$IUCN_ID,paste0(dir_niche,'/filtering_step3_IDs.rds'))
    
  }
  
  niches_f <- niches_f[niches_f$`Tmi_2.5%` > 273.15,]
  
  write.csv(niches_f,paste0(dir_niche,'niches_filtered.csv'),row.names = F)
  
  cat(paste0(rep('-',30)),'\n\n')
  
  
}

libinv(c('ggplot2','purrr','tidyr','bestNormalize','foreach'))

tab <- niches_f %>%
  as_tibble()

apply(tab %>% select(-IUCN_ID),2,function(x) sum(is.na(x)))
apply(tab %>% select(-IUCN_ID),2,function(x) sum(x == 0))


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

