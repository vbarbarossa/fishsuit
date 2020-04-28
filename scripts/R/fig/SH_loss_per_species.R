source('R/4targets/MASTER.R')

library(foreach); library(ggplot2); library(RColorBrewer)

ESH_filename <- paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds')

if(file.exists(ESH_filename)){
  tab <- readRDS(ESH_filename)
}else{
  
  tab <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    ESH_tab <- read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab.csv'))
    
    tcl <- reshape2::melt(ESH_tab,id.vars = c(1:2,ncol(ESH_tab)))
    colnames(tcl)[(ncol(tcl)-1):ncol(tcl)] <- c('ESH_type','ESH')
    levels(tcl$ESH_type) <- c('Total','Q','Tw','Q&Tw')
    tcl$GCM <- clmod
    return(tcl)
  }
  tab <- droplevels(tab)
  tab$GCM <- factor(tab$GCM)
  
  tab <- cbind(tab,as.data.frame(do.call('rbind',strsplit(as.character(tab[,3]),'_'))))
  colnames(tab)[7:9] <- c('RCP','warmt','year')
  
  saveRDS(tab,ESH_filename)
}

levels(tab$RCP) <- c('2.6','4.5','6.0','8.5')
levels(tab$GCM) <- c('GFDL','HadGEM','IPSL','MIROC','NorESM')

# select only total ESH
tab <- droplevels(tab[tab$ESH_type == 'Total',])
# for each species
# calculate average
t <- do.call('rbind',
             lapply(
               split(tab,tab$id_no),function(x){
                 do.call('rbind',
                         lapply(split(x,x$warmt),function(y){
                           if(nrow(y) > 0) return(data.frame(id_no = unique(y$id_no),SH_loss = mean(y$ESH,na.rm = T), warmt = unique(y$warmt)))
                         } 
                         )
                 )
               }
             ))

# make table 
library(tidyr)
t <- spread(t, warmt, SH_loss)

# merge with IUCN names
iucn <- rbind(foreign::read.dbf('data/FW_FISH_PART_1.dbf'),foreign::read.dbf('data/FW_FISH_PART_2.dbf')) %>%
  .[c('id_no','binomial')] %>% .[!duplicated(.),]
  
m <- droplevels(cbind(data.frame(binomial = merge(t,iucn,by = 'id_no')[,'binomial']),t[,2:5]))

write.csv(m,'fishsuit_completeRun_warming_4targets/tabs/SH_loss_per_species.csv') 


