#' ---
#' title: "Phylogenetic generalized least squares on species traits"
#' author: "Valerio Barbarossa"
#' date: "November 27, 2019"
#' output: html_document
#' ---
#' 
#' ## Data cleaning
#' 
#' Read the data (produced by stats_on_traits.R) and inspect
setwd('R/4targets')
valerioUtils::libinv('dplyr')
df <- read.csv('../../stats_on_traits_output/df_for_stats_on_traits.csv') %>%
  as_tibble()
df

#' read also the dispersal data output

df_dsp <- read.csv('../../stats_on_traits_output/df_for_stats_on_traits_dispersal.csv') %>%
  as_tibble()
df_dsp

#' and merge with df

# first select and rename the columns
df_dsp <- df_dsp %>%
  select(binomial,starts_with('RC'))
colnames(df_dsp) <- paste0(colnames(df_dsp),'_dsp')
df_dsp
# then merge
df <- left_join(df,df_dsp,by = c(binomial='binomial_dsp'))
df

#' 
#' check NAs
#' 
## ------------------------------------------------------------------------
apply(df,2,function(x) sum(is.na(x)))


#' 
#' * looks like foodtroph and importance have more than 5k NAs
#' * there are 34-8 NAs in the RC variables
#' * there are 524 NAs in length
#' 
#' Remove NA entries for length and RC and assign data deficient "DD" to NAs of importance and foodtrophcat
#' 
## ------------------------------------------------------------------------
df <- df %>%
  select(-FoodTroph,-order,-family) %>%
  filter(!is.na(RC1.5) & !is.na(length)) %>%
  mutate(foodtrophcat = as.character(foodtrophcat),
         importance = as.character(importance))

df$importance[is.na(df$importance)] <- 'DD'
df$foodtrophcat[is.na(df$foodtrophcat)] <- 'DD'

#' and check again

apply(df,2,function(x) sum(is.na(x)))

#' We want to estimate the rate of change for the range contractions.
#' 
#' Therefore for each species we fit a linear regression across the four points and **report its slope as the rate of change**.  
#' Litterally we are fitting RC = i + s*T + e for each species,
#' where RC is in %, i is the intercept, s is the slope, T is the warming target (1.5,2.0 ...) and e is the residual error

df$RC_rate <- apply(df[,paste0('RC',c('1.5','2.0','3.2','4.5'))],1,function(x) as.numeric(coef(lm(as.numeric(x) ~ c(1.5,2,3.2,4.5)))[2]))
df$RC_rate_dsp <- apply(df[,paste0('RC',c('1.5','2.0','3.2','4.5'),'_dsp')],1,function(x) as.numeric(coef(lm(as.numeric(x) ~ c(1.5,2,3.2,4.5)))[2]))

#' compare dispersal vs non dispersal rates
library(ggplot2)
ggplot(df) +
  geom_point(aes(x = RC_rate,y = RC_rate_dsp)) +
  geom_abline(slope=1,intercept=0,color = 'red')

#'
#' similar trend but a lot of variation, interesting!
#' 
#' ********************
#' 
#' ## Check variables
#' 
#' The **dataset** has `r nrow(df)` entries
#' 
#' Check distribution of continuos variables

valerioUtils::libinv(c('ggplot2','purrr','tidyr'))
df %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

#' area and length are very skewed -> log-transformation?
## ---- echo=FALSE---------------------------------------------------------
hist(log10(df$area))
hist(log10(df$length))

#' Much better! Then log-transform them in the dataset

df <- df %>%
  mutate(area = log10(area), length = log10(length))

#' check again distribution of variables

df %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

#' check correlation of covariates, with correlation matrix. 
#'    
#' Need to convert all categorical variables to binary first:
## ------------------------------------------------------------------------
# corrmatrix
cm <- cor(df %>%
            select(area,length, climate, habitat, code) %>% 
            fastDummies::dummy_cols(.) %>%
            select(-climate, -habitat, -code)
)

#' ..and then compute the correlation matrix

valerioUtils::libinv('corrplot')
corrplot(cm, method = 'number', type = 'lower',number.cex = 1)

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
#' Load required packages

valerioUtils::libinv(c(
  'car',
  'ape', # Paradis et al 2004
  'nlme', # regression modelling
  'phytools', # Revell et al 2012
  'geiger', #Harmon et al 2007
  'caper',
  'MASS',
  'stringr'
))

#' Load phylogenetic tree for bony fish species

fish.tre<-read.tree("../../betancurt-R et al 2017.tre")#reading phylogenetic tree (Phylogenetic classification of bony fishes)
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

#' This is actually pretty nice, shows the phylogeny for the `r nrow(fish2)` species available

plotTree(full.tree2,type="fan",fsize=0.2,lwd=1,ftype="i")

#' I am stopping with Felix example here.
#' 
#' **Dataset** now consists of `r nrow(fish2)` entries based on species available from the phylogenetic tree, I might need to look into synonyms later on..
#' 
#' Let's first of all test a simple bivariate correlation based on our data  
#' ..for RC3.2  
#' (sorry about the labels of the continuous variables!)

fish2 %>%
  as_tibble() %>%
  gather(-binomial, -species2, -starts_with('RC'), key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = RC3.2)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

#' ..and for RC_rate

fish2 %>%
  as_tibble() %>%
  gather(-binomial, -species2, -starts_with('RC'), key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = RC_rate)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

#' then we fit the PGLM model. There are two main of parameters to keep into account here
#' 
#' 1. corPagel is a variation of the brownian motion used to represent the phylogenentic evolution
#' 2. lambda parameter, which defines the error in the length of the tree branches length
#' + lambda = 0 (independence)
#' + lambda = 1 (Brownian motion)
#' 

## ------------------------------------------------------------------------
#' let the model self calibrate the lambda parameter

fitc <- gls(RC3.2 ~ area + length + habitat + climate + code + importance + foodtrophcat,
            correlation=corPagel(value=1, phy=full.tree2,fixed = FALSE),
            method = "ML",
            data=fish2)

summary(fitc)
Anova(fitc)

#' Lambda = `r round(as.numeric(fitc$modelStruct),2)`
#' 
#' let's check the residuals, although I am not sure we need to do it
#' 

plot(fitc, resid(., type="n")~fitted(.), main="Normalized Residuals v Fitted Values",
     abline=c(0,0))
res <- resid(fitc, type="n")
qqnorm(res)
qqline(res)


#' let's check the variable importance (special thanks to Mirza for providing the R fucntion that does it right away!)
#' for each RC variable

L <- as.numeric(fitc$modelStruct)

source('variable_importance.R') #function provided by Mirza
library(foreach)

for(var in colnames(df)[grep('RC',colnames(df))]){
  
  dfsel <- fish2
  dfsel$RC <- dfsel[,var]
  
  fit <- gls(RC ~ area + length + habitat + climate + code + importance + foodtrophcat,
             correlation=corPagel(value=L, phy=full.tree2,fixed = TRUE),
             method = "ML",
             data=dfsel)
  
  coefdf <- as.data.frame(summary(fit)$tTable) %>%
    mutate(variable = row.names(.)) %>%
    mutate(`p-value` = ifelse(`p-value` < 10**-3,'<0.001',round(`p-value`,3)),
           Value = round(Value,3)) %>%
    as_tibble() %>%
    dplyr::select(variable, coef = 'Value', pval = 'p-value')
  
  vi <- variable_importance(dfsel[,c('area','length','habitat','climate','code','importance','foodtrophcat')], fit, iterations_num = 100)
  
  vidf <- data.frame(
    var = row.names(vi),
    mean = apply(vi,1,function(x) mean(x)),
    sd = apply(vi,1,function(x) sd(x))
  ) %>%
    as_tibble()
  
  #rename the columns with the variable
  colnames(vidf)[2:3] <- paste0(var,'_',colnames(vidf)[2:3])
  colnames(coefdf)[2:3] <- paste0(var,'_',colnames(coefdf)[2:3])
  
  if(var == 'RC1.5'){
    tab_res <- vidf
    coef_res <- coefdf
    
  }else{
    tab_res = cbind(tab_res,vidf[,2:3])
    coef_res = cbind(coef_res,coefdf[,2:3])
  } 
  
}

tab_res

#' SDs of single response variables look negligible
#' 
#' Compute the mean across the RC of different warming targets and standard deviation for plotting
#' 

tab_res$RC_MEAN <- tab_res %>% dplyr::select(-contains('dsp'),-starts_with('RC_rate')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,mean)
tab_res$RC_SD <- tab_res %>% dplyr::select(-contains('dsp'),-starts_with('RC_rate')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,sd)
tab_res$RC_MEAN_dsp <- tab_res %>% dplyr::select(contains('dsp'),-starts_with('RC_rate')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,mean)
tab_res$RC_SD_dsp <- tab_res %>% dplyr::select(contains('dsp'),-starts_with('RC_rate')) %>% dplyr::select(contains('mean')) %>%
  apply(.,1,sd)

write.csv(tab_res,'../../fishsuit_completeRun_warming_4targets/figshare/phyloreg_var_importance.csv',row.names = F)

#' and plot it
#' 

# format the table for ggplot
dp <- foreach(v = as.character(tab_res$var),.combine = 'rbind') %do% {
  
  t <- tab_res %>%
    dplyr::select(var,starts_with('RC_MEAN'),starts_with('RC_SD'),starts_with('RC_rate')) %>%
    filter(var == v)
  
  data.frame(
    var = v,
    mean = c(t$RC_MEAN,t$RC_MEAN_dsp,t$RC_rate_mean,t$RC_rate_dsp_mean),
    sd = c(t$RC_SD,t$RC_SD_dsp,t$RC_rate_sd,t$RC_rate_dsp_sd),
    variable = c('Range Contractions ND','Range Contractions MD','Rate of change ND','Rate of change MD')
  )
  
} %>%
  as_tibble() %>%
  mutate(var = forcats::fct_recode(var,
                                   'Range area' = 'area',
                                   'Body length' = 'length',
                                   'Habitat type' = 'habitat',
                                   'Climate' = 'climate',
                                   'IUCN code' = 'code',
                                   'Commercial imp.' = 'importance',
                                   'Trophic cat.' = 'foodtrophcat'))


p <- ggplot(dp,aes(x=var, y=mean,fill=variable)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4,position = position_dodge(width=0.9)) +
  scale_fill_manual(values = viridis::viridis(10,option='C')[c(4,5,8,9)]) +
  xlab('') +
  ylab('Variable importance') +
  coord_flip(expand = F) +
  facet_grid(var~.,scales = 'free_y',space = 'free_y',switch = 'y') +
  guides(fill = guide_legend(title=NULL)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        legend.position = c(0.7,0.13),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1.2, "lines"),
        text = element_text(size = 15))

ggsave('../../fishsuit_completeRun_warming_4targets/figs/traits_barplots.jpg',p,width = 89,height = 100,units='mm',scale = 1.5,dpi = 600)
ggsave('../../fishsuit_completeRun_warming_4targets/figs/traits_barplots.pdf',p,width = 89,height = 100,units='mm',scale = 1.5)




