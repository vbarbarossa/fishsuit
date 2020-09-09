source('config.R'); # always load as functions are loaded within this script

# load packages
library(foreach); library(dplyr)

# make sure tabs directory exists
dir_('tabs')

for(suffix in c('','_dispersal2')){
  
  # compile master table with results from the modelling
  tab <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    t0 <- read.csv(paste0('proc/',clmod,'/ESH_tab',suffix,'.csv'))
    t <- cbind(t0,as.data.frame(do.call('rbind',strsplit(as.character(t0[,'comboscen']),'_'))))
    t$GCM <- clmod
    t$comboscen <- paste0(t$comboscen,'_',t$GCM)
    colnames(t)[colnames(t) %in% paste0('V',1:3)] <- c('RCP','warmt','year')
    
    return(t)
  }
  tab <- droplevels(tab)
  tab$GCM <- factor(tab$GCM)
  tab$comboscen <- factor(tab$comboscen)
  
  # for each scenario get the mean range loss for each flow and temperature variable
  tab_scen_mean <- lapply(split(tab %>% select(starts_with('ESH'),warmt,comboscen),tab$warmt),
                          function(x){
                            w <- x$warmt[1]
                            return(
                              x %>% group_by(comboscen) %>% 
                                summarise_all(function(x) mean(x,na.rm=T) ) %>%
                                mutate(warmt = w))
                          }) %>% do.call('rbind',.) %>% as_tibble()
  
  # average across warming targets
  tab_warmt_mean <- tab_scen_mean %>% select(-comboscen) %>% group_by(warmt) %>% 
    summarise_all(function(x){paste0(round(mean(x),2),' (',round(sd(x),2),')')})
  
  # and save
  write.csv(tab_warmt_mean, paste0('tabs/overall_range_loss_mean_by_warming_target',suffix,'.csv'),row.names = F)
  write.csv(tab_scen_mean, paste0('tabs/overall_range_loss_mean_by_scenario',suffix,'.csv'),row.names = F)
  
  # for each scenario get the portion of species losing more that 50% of range for each variable
  tab_scen_50 <- lapply(split(tab %>% select(starts_with('ESH'),warmt,comboscen),tab$warmt),
                         function(x){
                           w <- x$warmt[1]
                           return(
                             x %>% group_by(comboscen) %>% 
                               summarise_all(function(x){length(x[x >= 50 & !is.na(x)])/length(x) * 100}) %>%
                               mutate(warmt = w))
                         }) %>% do.call('rbind',.) %>% as_tibble()
  
  # average across warming targets
  tab_warmt_50 <- tab_scen_50 %>% select(-comboscen) %>% group_by(warmt) %>% 
    summarise_all(function(x){paste0(round(mean(x),2),' (',round(sd(x),2),')')})
  
  # and save
  write.csv(tab_warmt_50, paste0('tabs/overall_range_loss_maj50_by_warming_target',suffix,'.csv'),row.names = F)
  write.csv(tab_scen_50, paste0('tabs/overall_range_loss_maj50_by_scenario',suffix,'.csv'),row.names = F)
  
  
}
