source('config.R')

valerioUtils::libinv(c(
  'ape','stringr','geiger','dplyr','corrplot','ggplot2','foreach','nlme'
))

source('scripts/R/fun/variable_importance.R') #function provided by Mirza for variable importance

df <- read.csv('figshare/RC_by_species.csv') %>%
  as_tibble()
colnames(df) <- c('binomial',paste0('RC',rep(c('1.5','2.0','3.2','4.5'),each=2),c('_dsp','')))

df <- df %>%
  left_join(read.csv('proc/species_traits.csv')) %>%
  as_tibble()

# create one column for habitat
df$habitat <- 'Lo'
df$habitat[df$lotic == 1] <- 'Lo'
df$habitat[df$lentic == 1] <- 'Le'
df$habitat[df$marine == 1] <- 'Ma'

# adjust importance
levels(df$importance) <- c('Com.','Com.','Com.','NoInt.','NoInt.','Subs.')

# make climate zone a factor
df$climate_zone <- as.factor(df$climate_zone)
levels(df$climate_zone) <- LETTERS[1:5]

# make habitat a factor
df$habitat <- as.factor(df$habitat)

## ------------------------------------------------------------------------
df$code[is.na(df$code)] <- 'DD'
df$code[is.na(df$code)] <- 'DD' # 2,597 probably from new data that do not have a match in IUCN
df$code[df$code %in% c('EX','EW')] <- 'DD' # there are supposedly still 4 extinct species
df$code[df$code %in% c('LRlc')] <- 'LC'
df$code[df$code %in% c('LRnt')] <- 'NT'

df <- df %>%
  select(-FoodTroph,-lotic,-lentic,-lentic_only,-marine,-id_no) %>%
  filter(!is.na(length)) %>%
  mutate(foodtrophcat = as.character(foodtrophcat),
         importance = as.character(importance))

df$importance[is.na(df$importance)] <- 'DD'
df$foodtrophcat[is.na(df$foodtrophcat)] <- 'DD'

df$foodtrophcat <- as.factor(df$foodtrophcat)
df$importance <- as.factor(df$importance)

df <- droplevels(df)

cat('Dataset has:',nrow(df), 'species after filtering out some NAs in the traits')

df <- df %>%
  mutate(area = log10(area), length = log10(length))

fish <- as.data.frame(df)


ft <- readRDS('proc/fishtree_stochastic.rds')


#' Load phylogenetic tree for bony fish species
fish.tre <- ft[[1]]
# str(fish.tre)#Phylogenetic information for more than 11000 fish species!!!
tips<-fish.tre$tip.label#species labels contained in the tree

#' Adjust the tree
# is.binary.tree(fish.tre) # we want this to be TRUE, if FALSE, run line below
# is.ultrametric(fish.tre) # if not run line below
fish.tre<-compute.brlen(fish.tre,method = "Grafen")
# is.ultrametric(fish.tre) # check again
full.tree<-fish.tre

# checking list of species in the tree
fish$species2 <- str_replace_all(fish$binomial," " , '_')
rownames(fish)<-fish$species2 
# Do all the species in the data are in the tree & vice versa?.
nombres.full<-name.check(full.tree,fish)
full.tree2<-drop.tip(full.tree,nombres.full$tree_not_data)
names.full2<-name.check(full.tree2,fish)
exclude.sp<-as.data.frame(names.full2$data_not_tree)

fish2<-fish[ !(fish$species2 %in% exclude.sp$`names.full2$data_not_tree`), ]
#Now are all species in the tree contained in the dataframe and viceversa

# cat('Dataset has:',nrow(fish2), 'species with phylogenetic info')

qq_df <- data.frame(matrix(ncol=8,nrow = nrow(fish2)))
colnames(qq_df) <- colnames(df)[grep('RC',colnames(df))]


for(var in colnames(df)[grep('RC',colnames(df))]){
  
  fit_list <- lapply(list.files(paste0('proc/phyloreg_fit_stoch_',var,'/'),full.names = T),readRDS)
  
  # check CV of coefficients
  
  library("abind")
  coef_array <- abind(lapply(fit_list,function(x) summary(x)$tTable), along=3)
  
  ttab <- apply(coef_array, c(1,2), mean)
  cvtab <- apply(coef_array, c(1,2), function(x) sd(x)/mean(x)*100)
  signtab <- apply(coef_array, c(1,2), function(x) sum(sign(x))/length(fit_list))
  
  write.csv(cvtab,paste0(dir_('tabs/phyloreg_stochastic_diag/'),'cv_',var,'.csv'),row.names = F)
  write.csv(cvtab,paste0(dir_('tabs/phyloreg_stochastic_diag/'),'sign_',var,'.csv'),row.names = F)
  
  dfsel <- fish2
  dfsel$RC <- log10( ( dfsel[,var]*(nrow(dfsel)-1)+0.5 )/nrow(dfsel)  )
  
  # average residuals
  qq_df[,var] <- lapply(fit_list,function(x) residuals(x,type = 'pearson')) %>% Reduce('+',.)/length(fit_list)
  # qq_df[,var] <- residuals(fit_list[[1]],type = 'pearson')
  
  # average coefficients
  coefdf <- as.data.frame(ttab) %>%
    mutate(variable = row.names(.)) %>%
    mutate(`p-value` = ifelse(`p-value` < 10**-3,'<0.001',round(`p-value`,3)),
           Value = round(Value,3),
           tval = abs(`t-value`)) %>%
    as_tibble() %>%
    dplyr::select(variable, coef = 'Value', tval, pval = 'p-value')
  colnames(coefdf)[2:4] <- paste0(var,'_',colnames(coefdf)[2:4])
  
  # average VI
  vi_list <- lapply(fit_list, function(x){
    vi <- variable_importance(dfsel[,c('area','length','habitat','climate_zone','code','importance','foodtrophcat')], x, iterations_num = 100)
    vidf <- data.frame(
      var = row.names(vi),
      mean = apply(vi,1,function(x) mean(x)),
      sd = apply(vi,1,function(x) sd(x))
    )
    
    #rename the columns with the variable
    colnames(vidf)[2:3] <- paste0(var,'_',colnames(vidf)[2:3])
    return(vidf)
  } 
  )
  
  vi_array <- abind(lapply(vi_list,function(x) x[,-1]), along=3)
  
  vidf <- apply(vi_array, c(1,2), mean) %>% as_tibble() %>% mutate(var = vi_list[[1]]$var)
  visd <- apply(vi_array, c(1,2), sd)
  
  vidf[,grep('sd',colnames(vidf))] <- visd[,1]
  
  # record lambdas (useful when setting fixed = FALSE)
  lamb <- data.frame(a = c(as.numeric(fit_list[[1]]$modelStruct),cor(dfsel$RC,as.numeric(fit_list[[1]]$fitted))))
  colnames(lamb) <- var
  row.names(lamb) <- c('lambda','r')
  
  if(var == colnames(df)[grep('RC',colnames(df))][1]){
    tab_res <- cbind(vidf[,3],vidf[,1:2])
    coef_res <- coefdf
    lambdas <- lamb
    
  }else{
    tab_res = bind_cols(tab_res,vidf[,1:2])
    coef_res = bind_cols(coef_res,coefdf[,2:4])
    lambdas = bind_cols(lambdas,lamb)
  } 
  
}

tab_res
coef_res
lambdas

#' SDs of single response variables look negligible
#' 
#' Compute the mean across the RC of different warming targets and standard deviation for plotting
#' 

tab_res$RC_MEAN <- tab_res %>% dplyr::select(-contains('dsp')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,mean)
tab_res$RC_SD <- tab_res %>% dplyr::select(-contains('dsp')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,sd)
tab_res$RC_MEAN_dsp <- tab_res %>% dplyr::select(contains('dsp')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,mean)
tab_res$RC_SD_dsp <- tab_res %>% dplyr::select(contains('dsp')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,sd)

coef_res$RC_coef_MEAN <- coef_res %>% dplyr::select(-contains('dsp')) %>% dplyr::select(contains('coef')) %>%
  apply(.,1,mean)
coef_res$RC_coef_SD <- coef_res %>% dplyr::select(-contains('dsp')) %>% dplyr::select(contains('coef')) %>%
  apply(.,1,sd)
coef_res$RC_coef_MEAN_dsp <- coef_res %>% dplyr::select(contains('dsp')) %>% dplyr::select(contains('coef')) %>%
  apply(.,1,mean)
coef_res$RC_coef_SD_dsp <- coef_res %>% dplyr::select(contains('dsp')) %>% dplyr::select(contains('coef')) %>%
  apply(.,1,sd)

dir_('figshare'); dir_('tabs')
write.csv(tab_res,'figshare/phyloreg_var_importance_stochastic.csv',row.names = F)
write.csv(coef_res,'tabs/phyloreg_coefficients_stochastic.csv',row.names = F)
write.csv(lambdas,'tabs/phyloreg_lambdas_stochastic.csv',row.names = F)
write.csv(coef_res %>% select(-contains('dsp')),'tabs/phyloreg_coefficients_nodsp_stochastic.csv',row.names = F)
write.csv(coef_res %>% select(variable,contains('dsp')),'tabs/phyloreg_coefficients_dsp_stochastic.csv',row.names = F)

#' and plot it
#' 

# QQ PLOT OF STD RESIDUALS
# y = qnorm( c(0.25, 0.75))
# x = quantile(residuals(fit,type = 'pearson'), c(0.25,0.75), type = 5)
# 
# slope <- diff(y) / diff(x)
# int   <- y[1] - slope * x[1]
# qqnorm(fit, abline = c(int,slope), grid = T)
# 
(qqp <- ggplot(tidyr::gather(qq_df), aes(sample = value)) +
    stat_qq(size = 0.5) + stat_qq_line( color = 'blue') +
    ylab('Standardized residuals') +
    xlab('Quantiles of standard normal') +
    facet_wrap('key', ncol =2) +
    coord_cartesian(ylim = c(-4,4)) +
    theme_bw() +
    theme(strip.background = element_blank())
)

(qqp_nl <- ggplot(tidyr::gather(qq_df), aes(sample = value)) +
    stat_qq(size = 0.5) + 
    ylab('Standardized residuals') +
    xlab('Quantiles of standard normal') +
    facet_wrap('key', ncol =2) +
    coord_cartesian(ylim = c(-4,4)) +
    theme_bw() +
    theme(strip.background = element_blank())
)

ggsave('figs/phyloreg_qqplots_stochastic.jpg',qqp,width = 89,height = 89*2,units='mm',scale = 1,dpi = 600)
ggsave('figs/phyloreg_qqplots_noline_stochastic.jpg',qqp_nl,width = 89,height = 89*2,units='mm',scale = 1,dpi = 600)


# format the table for ggplot
dp <- foreach(v = as.character(tab_res$var),.combine = 'rbind') %do% {
  
  t <- tab_res %>%
    dplyr::select(var,starts_with('RC_MEAN'),starts_with('RC_SD')) %>%
    filter(var == v)
  
  data.frame(
    var = v,
    mean = c(t$RC_MEAN,t$RC_MEAN_dsp),
    sd = c(t$RC_SD,t$RC_SD_dsp),
    variable = c('No dispersal','Maximal dispersal')
  )
  
} %>%
  as_tibble() %>%
  mutate(var = forcats::fct_recode(var,
                                   'Range area' = 'area',
                                   'Body length' = 'length',
                                   'Habitat type' = 'habitat',
                                   'Climate' = 'climate_zone',
                                   'IUCN code' = 'code',
                                   'Commercial imp.' = 'importance',
                                   'Trophic cat.' = 'foodtrophcat'))


# insert lambda
levels(dp$variable) <- c(
  paste0('Max dispersal\nλ=',lambdas[1,] %>% select(contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-'),'; r=',lambdas[2,] %>% select(contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-')),
  paste0('No dispersal\nλ=',lambdas[1,] %>% select(-contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-'),'; r=',lambdas[2,] %>% select(-contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-'))
)

p <- ggplot(dp,aes(x=var, y=mean,fill=variable)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4,position = position_dodge(width=0.9)) +
  scale_fill_manual(values = viridis::viridis(10,option='C')[c(5,8)]) +
  xlab('') +
  ylab('Variable importance') +
  coord_flip(expand = F) +
  facet_grid(var~.,scales = 'free_y',space = 'free_y') +
  guides(fill = guide_legend(title=NULL,reverse = T)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        legend.position = c(0.7,0.15),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        text = element_text(size = 9))
p

ggsave('figs/traits_barplots_stochastic.jpg',p,width = 89,height = 80,units='mm',scale = 1,dpi = 600)
# ggsave('figs/traits_barplots.pdf',p,width = 89,height = 80,units='mm',scale = 1)

# do one per warming target
dp2 <- foreach(v = as.character(tab_res$var),.combine = 'rbind') %do% {
  
  t <- tab_res %>%
    dplyr::select(-starts_with('RC_MEAN'),-starts_with('RC_SD')) %>%
    filter(var == v)
  
  foreach(wt = warming_targets,.combine = 'rbind') %do%{
    
    data.frame(
      var = v,
      mean = c(as.data.frame(t)[,paste0('RC',wt,'_mean')],as.data.frame(t)[,paste0('RC',wt,'_dsp_mean')]),
      sd = c(as.data.frame(t)[,paste0('RC',wt,'_sd')],as.data.frame(t)[,paste0('RC',wt,'_dsp_sd')]),
      variable = c('No dispersal','Maximal dispersal'),
      wt = wt
    )
    
  }
  
} %>%
  as_tibble() %>%
  mutate(var = forcats::fct_recode(var,
                                   'Range area' = 'area',
                                   'Body length' = 'length',
                                   'Habitat type' = 'habitat',
                                   'Climate' = 'climate_zone',
                                   'IUCN code' = 'code',
                                   'Commercial imp.' = 'importance',
                                   'Trophic cat.' = 'foodtrophcat'))

levels(dp2$variable) <- c(
  paste0('Max dispersal\nλ=',lambdas[1,] %>% select(contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-'),'\nr=',lambdas[2,] %>% select(contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-')),
  paste0('No dispersal\nλ=',lambdas[1,] %>% select(-contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-'),'\nr=',lambdas[2,] %>% select(-contains('dsp')) %>% range %>% round(2) %>% paste0(collapse = '-'))
)

levels(dp2$wt) <- c('1.5°C','2.0°C','3.2°C','4.5°C')

p <- ggplot(dp2,aes(x=var, y=mean,fill=variable)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.5,position = position_dodge(width=0.9)) +
  scale_fill_manual(values = viridis::viridis(10,option='C')[c(5,8)]) +
  xlab('') +
  ylab('Variable importance') +
  coord_flip(expand = F) +
  facet_grid(var~wt,scales = 'free_y',space = 'free_y') +
  guides(fill = guide_legend(title=NULL,reverse = T)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        legend.position = c(0.92,0.2),
        strip.text.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed',color='black'),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        text = element_text(size = 12))
p

ggsave('figs/traits_barplots_warming_levels_stochastic.jpg',p,width = 89*4,height = 80*2,units='mm',scale = 0.6,dpi = 600)

ggsave('figs/traits_barplots_warming_levels_stochastic_l.jpg',p,width = 213.6,height = 96,units='mm',dpi = 600)
ggsave('figs/traits_barplots_warming_levels_stochastic_l.pdf',p,width = 210,height = 95,units='mm',dpi = 600)

ggsave('figs/traits_barplots_warming_levels_stochastic.jpg',p,width = 89*4,height = 80*2,units='mm',scale = 0.6,dpi = 600)
# ggsave('figs/traits_barplots.pdf',p,width = 89,height = 80,units='mm',scale = 1)

