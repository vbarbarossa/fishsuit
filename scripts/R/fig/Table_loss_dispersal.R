source('R/4targets/MASTER.R')

library(foreach); library(ggplot2); library(RColorBrewer)

ESH_filename <- 'tabs/ESH_merged_warmingtargets_confidencemargin_dispersal.rds'

if(file.exists(ESH_filename)){
  tab <- readRDS(ESH_filename)
}else{
  
  tab <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    t0 <- read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab_dispersal.csv'))
    t <- cbind(t0,as.data.frame(do.call('rbind',strsplit(as.character(t0[,'comboscen']),'_'))))
    t$GCM <- clmod
    t$comboscen <- paste0(t$comboscen,'_',t$GCM)
    colnames(t)[8:10] <- c('RCP','warmt','year')
    
    return(t)
  }
  tab <- droplevels(tab)
  tab$GCM <- factor(tab$GCM)
  tab$comboscen <- factor(tab$comboscen)
  saveRDS(tab,ESH_filename)
}

wt_tab <- list()
c=1
summary_tab <- foreach(wt = levels(tab$warmt),.combine = 'cbind') %do% {
  
  tt <- droplevels(tab[tab$warmt == wt,])
  
  sc_tab <- foreach(sc = levels(tt$comboscen),.combine = 'cbind') %do% {
    t <- droplevels(tt[tt$comboscen == sc,])
    v <- t[!is.na(t$ESH_all),'ESH_all']
    v <- v[!is.infinite(v)]
    c(
      length(v),median(v),
      length(v[v >= 0])/length(v)*100,
      length(v[v >= 25])/length(v)*100,
      length(v[v >= 50])/length(v)*100,
      length(v[v < 0])/length(v)*100,
      length(v[v < -50])/length(v)*100,
      length(v[v < -100])/length(v)*100,
      length(v[v < -500])/length(v)*100,
      length(v[v < -1000])/length(v)*100,
      length(v[v < -5000])/length(v)*100
      )
  }
  colnames(sc_tab) <- levels(tt$comboscen)
  row.names(sc_tab) <- c(
    'N','median',paste0('%>=',c(0,25,50)),paste0('%<',c(0,-50,-100,-500,-1000,-5000))
  )
  wt_tab[[c]] <- as.data.frame(apply(sc_tab,1,mean))
  c=c+1
  return(sc_tab)
}
wt_tab <- do.call('cbind',wt_tab)
colnames(wt_tab) <- warming_targets

write.csv(summary_tab,paste0(dir_mod,'tabs/SHL_summary_stats_per_scen_dispersal.csv'))
write.csv(wt_tab,paste0(dir_mod,'tabs/SHL_summary_stats_per_wt_dispersal.csv'))

# present-day dispersal
wt_tab <- list()
c=1
summary_tab <- foreach(wt = levels(tab$warmt),.combine = 'cbind') %do% {
  
  tt <- droplevels(tab[tab$warmt == wt,])
  
  sc_tab <- foreach(sc = levels(tt$comboscen),.combine = 'cbind') %do% {
    t <- droplevels(tt[tt$comboscen == sc,])
    v <- t[!is.na(t$ESH_all),'ESH_all_hist']
    v <- v[!is.infinite(v)]
    c(
      length(v),median(v),
      length(v[v >= 0])/length(v)*100,
      length(v[v >= 25])/length(v)*100,
      length(v[v >= 50])/length(v)*100,
      length(v[v < 0])/length(v)*100,
      length(v[v < -50])/length(v)*100,
      length(v[v < -100])/length(v)*100,
      length(v[v < -500])/length(v)*100,
      length(v[v < -1000])/length(v)*100,
      length(v[v < -5000])/length(v)*100
    )
  }
  colnames(sc_tab) <- levels(tt$comboscen)
  row.names(sc_tab) <- c(
    'N','median',paste0('%>=',c(0,25,50)),paste0('%<',c(0,-50,-100,-500,-1000,-5000))
  )
  wt_tab[[c]] <- as.data.frame(apply(sc_tab,1,mean))
  c=c+1
  return(sc_tab)
}
wt_tab <- do.call('cbind',wt_tab)
colnames(wt_tab) <- warming_targets

write.csv(summary_tab,paste0(dir_mod,'tabs/SHL_summary_stats_per_scen_dispersal_hist.csv'))
write.csv(wt_tab,paste0(dir_mod,'tabs/SHL_summary_stats_per_wt_dispersal_hist.csv'))


#_---------------------------------------------------------------------------------
# lotic only

ESH_filename <- 'tabs/ESH_merged_warmingtargets_confidencemargin_lotic.rds'

if(file.exists(ESH_filename)){
  tab <- readRDS(ESH_filename)
}else{
  
  tab <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    t0 <- read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab_lotic.csv'))
    t <- cbind(t0,as.data.frame(do.call('rbind',strsplit(as.character(t0[,'comboscen']),'_'))))
    t$GCM <- clmod
    t$comboscen <- paste0(t$comboscen,'_',t$GCM)
    colnames(t)[8:10] <- c('RCP','warmt','year')
    
    return(t)
  }
  tab <- droplevels(tab)
  tab$GCM <- factor(tab$GCM)
  tab$comboscen <- factor(tab$comboscen)
  saveRDS(tab,ESH_filename)
}


wt_tab <- list()
c=1
summary_tab <- foreach(wt = levels(tab$warmt),.combine = 'cbind') %do% {
  
  tt <- droplevels(tab[tab$warmt == wt,])
  
  sc_tab <- foreach(sc = levels(tt$comboscen),.combine = 'cbind') %do% {
    t <- droplevels(tt[tt$comboscen == sc,])
    v <- t[!is.na(t$ESH_all),'ESH_all']
    v <- v[!is.infinite(v)]
    c(
      length(v),median(v),
      length(v[v >= 0])/length(v)*100,
      length(v[v >= 25])/length(v)*100,
      length(v[v >= 50])/length(v)*100,
      length(v[v < 0])/length(v)*100,
      length(v[v < -50])/length(v)*100,
      length(v[v < -100])/length(v)*100,
      length(v[v < -500])/length(v)*100,
      length(v[v < -1000])/length(v)*100,
      length(v[v < -5000])/length(v)*100
    )
  }
  colnames(sc_tab) <- levels(tt$comboscen)
  row.names(sc_tab) <- c(
    'N','median',paste0('%>=',c(0,25,50)),paste0('%<',c(0,-50,-100,-500,-1000,-5000))
  )
  wt_tab[[c]] <- as.data.frame(apply(sc_tab,1,mean))
  c=c+1
  return(sc_tab)
}
wt_tab <- do.call('cbind',wt_tab)
colnames(wt_tab) <- warming_targets

write.csv(summary_tab,paste0(dir_mod,'tabs/SHL_summary_stats_per_scen_lotic.csv'))
write.csv(wt_tab,paste0(dir_mod,'tabs/SHL_summary_stats_per_wt_lotic.csv'))

#_---------------------------------------------------------------------------------
