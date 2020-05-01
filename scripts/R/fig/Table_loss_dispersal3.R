source('R/4targets/MASTER.R')

source('R/functions.R'); source('R/_dir_loc.R')

library(foreach); library(ggplot2); library(RColorBrewer)

climate_models <- c('gfdl','ipsl','hadgem','miroc','noresm')

ESH_filename <- 'tabs/ESH_merged_warmingtargets_confidencemargin_dispersal.rds'

if(file.exists(ESH_filename)){
  tab <- readRDS(ESH_filename)
}else{
  
  tab <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    t0 <- read.csv(paste0('proc/',clmod,'/ESH_tab_dispersal.csv'))
    t <- cbind(t0,as.data.frame(do.call('rbind',strsplit(as.character(t0[,'comboscen']),'_'))))
    t$GCM <- clmod
    t$comboscen <- paste0(t$comboscen,'_',t$GCM)
    colnames(t)[which(colnames(t) == 'V1'):which(colnames(t) == 'V3')] <- c('RCP','warmt','year')
    
    return(t)
  }
  tab <- droplevels(tab)
  tab$GCM <- factor(tab$GCM)
  tab$comboscen <- factor(tab$comboscen)
  saveRDS(tab,ESH_filename)
}

summary_tab <- foreach(wt = levels(tab$warmt),.combine = 'cbind') %do% {
  
  tt <- droplevels(tab[tab$warmt == wt,])
  
  sc_tab <- foreach(sc = levels(tt$comboscen),.combine = 'cbind') %do% {
    t <- droplevels(tt[tt$comboscen == sc,])
    v <- t[!is.na(t$ESH_all),'ESH_all']
    c(
      length(v),mean(v),median(v),IQR(v),
      length(v[v == 0])/length(v)*100,
      length(v[v >= 5])/length(v)*100,
      length(v[v >= 25])/length(v)*100,
      length(v[v >= 50])/length(v)*100,
      length(v[v >= 75])/length(v)*100,
      length(v[v >= 95])/length(v)*100,
      length(v[v >= 100])/length(v)*100,
      quantile(v,c(seq(0.025,0.975,0.025)))
    )
  }
  colnames(sc_tab) <- levels(tt$comboscen)
  row.names(sc_tab) <- c(
    'N','mean','median','IQR','%=0',paste0('%>=',c(5,25,50,75,95,100)),paste0('q',seq(0.025,0.975,0.025)*100)
  )
  return(sc_tab)
}

write.csv(summary_tab,paste0(dir_mod,'tabs/SHL_summary_stats_per_scen_dispersal.csv'))

#_---------------------------------------------------------------------------------

# this is available from 'Figure_SI_boxplot_scenario.R'
tab <- readRDS(paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds'))
tab <- droplevels( tab[tab$ESH_type == 'Total',] )

# custom_targets <- c('1.5','2.0','3.0')


tab_median <- foreach(sp = unique(tab$id_no),.combine = 'rbind') %do% {
  tsp <- droplevels(tab[tab$id_no == sp,])
  foreach(wt = levels(tab$warmt),.combine = 'rbind') %do% {
    
    t <- droplevels(tsp[tsp$warmt == wt,])
    return(
      data.frame(
        id_no = sp,
        no.cells = unique(t$no.cells),
        ESH_median = median(t$ESH,na.rm=T),
        warmt = wt
      )
    )
  } 
}


summary_tab <- foreach(wt = levels(tab_median$warmt),.combine = 'cbind') %do% {
  
  t <- droplevels(tab_median[tab_median$warmt == wt,])
  v <- t[!is.na(t$ESH_median),'ESH_median']
  c(
    length(v),mean(v),median(v),IQR(v),
    length(v[v == 0])/length(v)*100,
    length(v[v >= 5])/length(v)*100,
    length(v[v >= 25])/length(v)*100,
    length(v[v >= 50])/length(v)*100,
    length(v[v >= 75])/length(v)*100,
    length(v[v >= 95])/length(v)*100,
    length(v[v >= 100])/length(v)*100,
    quantile(v,c(seq(0.025,0.975,0.025)))
  )
}
colnames(summary_tab) <- levels(tab_median$warmt)
row.names(summary_tab) <- c(
  'N','mean','median','IQR','%=0',paste0('%>=',c(5,25,50,75,95,100)),paste0('q',seq(0.025,0.975,0.025)*100)
)

write.csv(summary_tab,paste0(dir_mod,'tabs/SHL_summary_stats.csv'))
