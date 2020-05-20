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

#' 
#' check NAs
#' 
## ------------------------------------------------------------------------
apply(df,2,function(x) sum(is.na(x)))


#' 
#' * looks like foodtroph and importance have more than 5k NAs
#' 
#' * there are 695 NAs in length
#' 
#' * there are 3522 NAs in code
#' 
#' Remove NA entries for length and RC and assign data deficient "DD" to NAs of importance and foodtrophcat
#' 
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

#' and check again

apply(df,2,function(x) sum(is.na(x)))

cat('Dataset has:',nrow(df), 'species after filtering out some NAs in the traits')

#' Much better! Then log-transform them in the dataset
df <- df %>%
  mutate(area = log10(area), length = log10(length))

## ------------------------------------------------------------------------
# corrmatrix
cm1 <- cor(df %>%
            select(area,length, climate_zone, habitat, code, importance, foodtrophcat) %>% 
            fastDummies::dummy_cols(.) %>%
            select(-climate_zone, -habitat, -code, -importance, -foodtrophcat)
)

cm2 <- cor(df %>% 
            select('Range area' = area,'Body length' = length, 'Climate' = climate_zone, 
                    'Habitat type' = habitat, 'IUCN code' = code, 'Commercial imp.' = importance, 
                   'Trophic cat.' = foodtrophcat) %>% 
            mutate_all(as.numeric), method = 'spearman'
)

# library(psych)
# cm <- mixedCor(df %>% 
#                  select(area,length, climate_zone, habitat, code, importance, foodtrophcat) %>% 
#                  mutate_all(as.numeric), 
#                method = 'pearson')

#' ..and then compute the correlation matrix
dir_('figs')
jpeg('figs/phyloreg_CORRPLOT.jpg',width = 300,height = 300,units='mm',res = 600,type='cairo')
corrplot(cm1, method = 'number', type = 'lower',number.cex = 1)
dev.off()

jpeg('figs/phyloreg_CORRPLOT_discretized.jpg',width = 150,height = 150,units='mm',res = 600,type='cairo')
corrplot(cm2, method = 'number', type = 'lower',number.cex = 1)
dev.off()

#' and with variance inflation factors (VIFs)
#' 
## ------------------------------------------------------------------------

# valerioUtils::libinv('HH')
# vif(as.data.frame(df %>%
#                     dplyr::select(area,length, climate, habitat, code) %>% 
#                     fastDummies::dummy_cols(.) %>%
#                     dplyr::select(-climate, -habitat, -code)
# ))

#' apparently cannot apply VIF to binary variables or factors, but the correlation plot was convincing enough!
#'  
#' ********************
#' 
#' ## Phylogenetic generalized linear models
#' 
#' Convert tibble to data.frame (apparently needed to use the PGLMs packages)
#' 
## ------------------------------------------------------------------------
fish <- as.data.frame(df)

#' Based on example suggested by Felix Leiva (RU):
#' 
#' Load phylogenetic tree for bony fish species

fish.tre<-read.tree("data/betancurt-R et al 2017.tre")#reading phylogenetic tree (Phylogenetic classification of bony fishes)
str(fish.tre)#Phylogenetic information for more than 11000 fish species!!!
tips<-fish.tre$tip.label#species labels contained in the tree
fish.tre$Nnode#Node numbers

#' Adjust the tree

is.binary.tree(fish.tre) # we want this to be TRUE, if FALSE, run line below
is.ultrametric(fish.tre) # if not run line below
fish.tre<-compute.brlen(fish.tre,method = "Grafen")
is.ultrametric(fish.tre) # check again
full.tree<-fish.tre

#' Quite redundant piece of code to exclude species not in the phylotree, since it works I am keeping it as it is for now

#checking list of species in the tree
fish$species2 <- str_replace_all(fish$binomial," " , '_')
rownames(fish)<-fish$species2 

# Do all the species in the data are in the tree & vice versa?.
nombres.full<-name.check(full.tree,fish)


full.tree2<-drop.tip(full.tree,nombres.full$tree_not_data)
# full.tree2$tip.label# 2,757 species in the data base with phylogenetic information

names.full2<-name.check(full.tree2,fish)
# But still there are 3,614 species contained in the dataframe but NO in the tree
exclude.sp<-as.data.frame(names.full2$data_not_tree)

fish2<-fish[ !(fish$species2 %in% exclude.sp$`names.full2$data_not_tree`), ]
#Now are all species in the tree contained in the dataframe and viceversa

#' then we fit the PGLM model. There are two main of parameters to keep into account here
#' 
#' 1. corPagel is a variation of the brownian motion used to represent the phylogenentic evolution
#' 2. lambda parameter, which defines the error in the length of the tree branches length
#' + lambda = 0 (independence)
#' + lambda = 1 (Brownian motion)
#' 

cat('Dataset has:',nrow(fish2), 'species with phylogenetic info')

# define settings for Pagel's covariance matrix
CP <- corPagel(value=1, phy=full.tree2,fixed = FALSE)

for(var in colnames(df)[grep('RC',colnames(df))]){
  
  dfsel <- fish2
  dfsel$RC <- dfsel[,var]
  
  fit <- gls(RC ~ area + length + habitat + climate_zone + code + importance + foodtrophcat,
             correlation=CP,
             method = "ML",
             data=dfsel)
  
  coefdf <- as.data.frame(summary(fit)$tTable) %>%
    mutate(variable = row.names(.)) %>%
    mutate(`p-value` = ifelse(`p-value` < 10**-3,'<0.001',round(`p-value`,3)),
           Value = round(Value,3),
           tval = abs(`t-value`)) %>%
    as_tibble() %>%
    dplyr::select(variable, coef = 'Value', tval, pval = 'p-value')
  
  vi <- variable_importance(dfsel[,c('area','length','habitat','climate_zone','code','importance','foodtrophcat')], fit, iterations_num = 100)
  
  vidf <- data.frame(
    var = row.names(vi),
    mean = apply(vi,1,function(x) mean(x)),
    sd = apply(vi,1,function(x) sd(x))
  ) %>%
    as_tibble()
  
  #rename the columns with the variable
  colnames(vidf)[2:3] <- paste0(var,'_',colnames(vidf)[2:3])
  colnames(coefdf)[2:4] <- paste0(var,'_',colnames(coefdf)[2:4])
  
  # record lambdas (useful when setting fixed = FALSE)
  lamb <- data.frame(a = c(as.numeric(fit$modelStruct),cor(dfsel$RC,as.numeric(fit$fitted))))
  colnames(lamb) <- var
  row.names(lamb) <- c('lambda','r')
  
  if(var == colnames(df)[grep('RC',colnames(df))][1]){
    tab_res <- vidf
    coef_res <- coefdf
    lambdas <- lamb
    
  }else{
    tab_res = bind_cols(tab_res,vidf[,2:3])
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
write.csv(tab_res,'figshare/phyloreg_var_importance.csv',row.names = F)
write.csv(coef_res,'tabs/phyloreg_coefficients.csv',row.names = F)
write.csv(lambdas,'tabs/phyloreg_lambdas.csv',row.names = F)
write.csv(coef_res %>% select(-contains('dsp')),'tabs/phyloreg_coefficients_nodsp.csv',row.names = F)
write.csv(coef_res %>% select(variable,contains('dsp')),'tabs/phyloreg_coefficients_dsp.csv',row.names = F)

#' and plot it
#' 

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

ggsave('figs/traits_barplots.jpg',p,width = 89,height = 80,units='mm',scale = 1,dpi = 600)
# ggsave('figs/traits_barplots.pdf',p,width = 89,height = 80,units='mm',scale = 1)

