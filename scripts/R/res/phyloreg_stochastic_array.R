g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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

# produce the fishtree list
# ft <- fishtree::fishtree_complete_phylogeny(fish$binomial)
# saveRDS(ft,'proc/fishtree_stochastic.rds')

var <- colnames(df)[grep('RC',colnames(df))][g]
ft <- readRDS('proc/fishtree_stochastic.rds')
dir_save <- dir_(paste0('proc/phyloreg_fit_stoch_',var,'/'))

fit_model <- function(j){
  
  cat('\nRepetition ',j,'\n')
  
  #' Load phylogenetic tree for bony fish species
  fish.tre <- ft[[j]]
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
  
  dfsel <- fish2
  # response variable transformation (due to skewness in qq plot)
  # 1. Smithson and Verkuilen to remove zeroes y'=(y*(n-1)+0.5)/n 
  # 2. log10
  dfsel$RC <- log10( ( dfsel[,var]*(nrow(dfsel)-1)+0.5 )/nrow(dfsel)  )
 
  # define settings for Pagel's covariance matrix
  if(j == 1){
    CP <- corPagel(value=1, phy=full.tree2,fixed = FALSE)
  }else{
    fit1 <- readRDS(paste0('proc/phyloreg_fit_stoch_',var,'/1.rds'))
    lamb_val <- as.numeric(fit1$modelStruct)
    CP <- corPagel(value=lamb_val, phy=full.tree2,fixed = TRUE)
  }
  
  fit <- gls(RC ~ area + length + habitat + climate_zone + code + importance + foodtrophcat,
             correlation=CP,
             method = "ML",
             data=dfsel)
  
  saveRDS(fit,paste0(dir_save,j,'.rds'))
  
}

# first calibrate lambda based on the first repetition
# fit_model(1) # <<<<<<<<<<momentarily commented as the first run had gone through, the rest not..

# then run in parallel the rest
parallel::mcmapply(fit_model,as.list(2:100),SIMPLIFY = FALSE,mc.cores = 5)


