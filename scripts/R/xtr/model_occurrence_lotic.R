#sbatch --array=1-5

slurm_arrayid<-Sys.getenv("SLURM_ARRAY_TASK_ID")
nodenr<-as.numeric(slurm_arrayid)

for(g in nodenr:nodenr){
  
  source('config.R'); # always load as functions are loaded within this script
  
  filter_lentic_out <- TRUE #override config.R
  filter_lotic_out <- FALSE
  
  libinv(c('raster','foreach'))
  
  clmod = climate_models[g]
  scen = scenarios[scenarios != 'hist'] #<
  years = floor(mean(timespan_scen)) #<
  
  
  cat(paste0(rep('-',30)),'\n')
  cat('thresholds used:\n')
  cat(paste(vars,parse=' '))
  cat('\n')
  cat(paste(thresholds,parse=' '))
  cat('\n',paste0(rep('-',30)),'\n\n')
  
  # # # # load iucn data
  iucn <- foreach(i = 1:2,.combine='rbind') %do% foreign::read.dbf(paste0(dir_data,'/FW_FISH_20181113/FW_FISH_PART_',i,'.dbf'))
  iucn <- iucn[,c('id_no','binomial')][!duplicated(iucn[,c('id_no','binomial')]),]
  row.names(iucn) <- 1:nrow(iucn)
  
  # # # # load species habitat type <<<<<<< need to check why we are losing almost 100 species when merging
  hab <- merge(iucn,read.csv(paste0(dir_master,'iucn_habitat_type.csv')),by='binomial')
  
  niche_list <- list()
  for(j in seq_along(vars)){
    niche <- readRDS(paste0(dir_model,clmod,'/niches/',vars[j],'.rds'))
    
    # remove ranges with NAs
    # remove ranges with SD = 0, as they are sampled for only one 5 arcmin cell
    # niche <- niche[!is.na(niche$mean),] # & niche$sd != 0
    niche <- niche[,c('IUCN_ID',paste0(thresholds[j],'%'))]
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
  # saveRDS(niche$IUCN_ID,paste0(dir_model,clmod,'/niches/filtering_step0_IDs.rds'))
  
  # STEP 1 ###
  # filter out species outside the pcrglobwb domain and
  # with less than a user defined number of grid-cells threshold
  d <- readRDS(paste0(dir_model,clmod,'/niches/Qmi.rds'))
  step1 <- d[!is.na(d$mean) & d$no.grids >= min_no_grid_cells,1]
  
  cat('STEP 1: removing ',(nrow(niche) - length(step1)),' species represented by less than ',min_no_grid_cells,' grid-cell or falling out of the pcr-globwb domain\n')
  niche <- niche[niche$IUCN_ID %in% step1,]
  # saveRDS(niche$IUCN_ID,paste0(dir_model,clmod,'/niches/filtering_step1_IDs.rds'))
  
  niches_f <- niche
  
  if(filter_lentic_out){
    
    # STEP 2
    ## merge habitat
    niches <- merge(niche,hab[,c('id_no','binomial','lotic','lentic')],by.x='IUCN_ID',by.y='id_no')
    cat('removing ',(nrow(niche) - nrow(niches)),' species when merging with habitat info\n')
    # saveRDS(niches$IUCN_ID,paste0(dir_model,clmod,'/niches/filtering_step2_IDs.rds'))
    
    # STEP 3
    ## filter for lotic sp only
    niches_f <- niches[(niches$lotic == 1),]
    cat('removing ',(nrow(niches) - nrow(niches_f)),' species exclusively lentic\n')
    # saveRDS(niches_f$IUCN_ID,paste0(dir_model,clmod,'/niches/filtering_step3_IDs.rds'))
    write.csv(niches_f,paste0(dir_model,clmod,'/niches/niches_filtered_lotic.csv'))
  }
  
  if(filter_lotic_out){
    
    # STEP 2
    ## merge habitat
    niches <- merge(niche,hab[,c('id_no','binomial','lotic','lentic')],by.x='IUCN_ID',by.y='id_no')
    cat('removing ',(nrow(niche) - nrow(niches)),' species when merging with habitat info\n')
    # saveRDS(niches$IUCN_ID,paste0(dir_model,clmod,'/niches/filtering_step2_IDs.rds'))
    
    # STEP 3
    ## filter for lotic sp only
    niches_f <- niches[(niches$lotic == 0),]
    cat('removing ',(nrow(niches) - nrow(niches_f)),' lotic species\n')
    # saveRDS(niches_f$IUCN_ID,paste0(dir_model,clmod,'/niches/filtering_step3_IDs.rds'))
    write.csv(niches_f,paste0(dir_model,clmod,'/niches/niches_filtered_lentic.csv'))
  }
  
  
  
  cat(paste0(rep('-',30)),'\n\n')
  
  #-------------------------------------------------------------------------------------------------------
  
  # make tab niches
  nich <- droplevels(niches_f[,1:(length(vars)+1)])
  colnames(nich) <- c('id_no',paste0(vars,'_th'))
  
  dir_out <- dir_(paste0(dir_model,clmod,'/modelled_occurrence//'))
  
  attr.tab <- foreach(i = seq_along(warming_targets),.combine = 'rbind') %do% {
    data.frame(
      clmod = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[1]),
      scen = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[2]),
      warmt = paste0(warming_targets[i],'C'),
      year = warming_tab[,i+1]
    )
  }
  
  attr.tab <- attr.tab[!is.na(attr.tab$year) & attr.tab$clmod == clmod,]
  row.names(attr.tab) <- NULL
  
  comboscen <- as.character(apply(attr.tab[,2:ncol(attr.tab)],1,function(x) paste(x,collapse = '_')))
  
  # ESH_tab
  cat(paste0(rep('-',30)),'\n')
  cat('Build ESH_tab\n')
  
  tabs <- list.files(dir_out,pattern='.rds')
  ids <- as.character(nich$id_no)
  
  
  # function to calculate ESH for each species-tab
  esh_calc <- function(sp){ #,sn
    # read the table
    t <- readRDS(paste0(dir_out,sp,'.rds'))
    t <- t[complete.cases(t),] #remove rows that contain one or more NAs
    
    d <- data.frame(
      id_no = rep(sp,length(comboscen)),
      no.cells = nrow(t),
      ESH_all = NA,
      ESH_Q = NA,
      ESH_T = NA,
      ESH_QnT = NA,
      comboscen = comboscen #sn is the looped scen
    )
    for(y in 1:nrow(d)){
      
      sn = d$comboscen[y]
      
      # select scenario and year
      th <- t[,paste0(vars,'_hist')] #historical
      ts <- t[,paste0(vars,'_',sn)] #scenario
      #filter out rows of ts that are already excluded in th due to the quantile threshold
      ts <- ts[-which(unname(apply(th,1,all) == FALSE)),]
      
      d$ESH_all[y] <- sum(apply(ts,1,all) == FALSE)/nrow(ts)*100
      d$ESH_Q[y] <- sum(apply(ts[,paste0(c('Qmi','Qzf'),'_',sn)],1,all) == FALSE)/nrow(ts)*100
      d$ESH_T[y] <- sum(ts[,paste0(c('Tma'),'_',sn)] == FALSE)/nrow(ts)*100
      d$ESH_QnT[y] <- sum(apply(cbind(apply(ts[,paste0(c('Qmi','Qzf'),'_',sn)],1,all),
                                      ts[,paste0(c('Tma'),'_',sn)])
                                ,1,any) == FALSE)/nrow(ts)*100
    }
    return(d)
  }
  
  # calculate ESH for each scenario and year span
  ESH_tab <- 
    do.call('rbind',
            parallel::mcmapply(esh_calc,sp = ids,SIMPLIFY = F,mc.cores = ncores)
    )
  
  row.names(ESH_tab) <- NULL
  ESH_tab$comboscen <- as.factor(ESH_tab$comboscen)
  
  if(filter_lentic_out) out <- paste0(dir_model,clmod,'/ESH_tab_lotic.csv')
  if(filter_lotic_out) out <- paste0(dir_model,clmod,'/ESH_tab_lentic.csv')
  
  write.csv(ESH_tab,out,row.names = F)
  
  cat('successfully written to ',out,'\n')
  cat(paste0(rep('-',30)),'\n\n')
  
  
  #----------------------------------------------------------------------------------------------------
  # SUMMARIZE at the cell-level
  
  cat(paste0(rep('-',30)),'\n')
  cat('Build SR_tab\n')
  
  parallel::mcmapply(function(sc) {
    
    #<<<<<< need to remember how this was done and re-do it (maybe in the initial scripts 01-03 on milkun )
    tab <- readRDS(paste0(dir_master,'ssp/tab_template.rds'))
    
    for(v in c(vars,'Q_all','T_all','any','all','both_QT')) tab@data[,v] <- NA
    
    # <<<<<<<< need to parallelize this
    for(i in ids){
      
      t <- readRDS(paste0(dir_out,i,'.rds'))
      t$X <- as.numeric(row.names(t)) #the way it was coded based on csv
      th <- t[,paste0(vars,'_hist')] #historical
      
      t <- t[-which(unname(apply(th,1,all) == FALSE)),]
      
      #all
      tab@data$occ[t$X] <- apply(cbind(tab@data$occ[t$X],rep(1,nrow(t))),1,function(x) sum(x,na.rm=T)) 
      
      tab@data[t$X,'all'] <- apply(cbind(tab@data[t$X,'all'],as.integer(apply(t[,paste0(c('Qmi','Qzf','Tma'),'_',sc)],1,all)))
                                   ,1,function(x) sum(x,na.rm=T))
      
      
      tab@data[t$X,'any'] <- apply(cbind(tab@data[t$X,'any'],as.integer(apply(t[,paste0(c('Qmi','Qzf','Tma'),'_',sc)],1,any)))
                                   ,1,function(x) sum(x,na.rm=T))
      
      
      tab@data[t$X,'Q_all'] <- apply(cbind(tab@data[t$X,'Q_all'],as.integer(apply(t[,paste0(c('Qmi','Qzf'),'_',sc)],1,all)))
                                     ,1,function(x) sum(x,na.rm=T))
      
      
      tab@data[t$X,'T_all'] <- apply(cbind(tab@data[t$X,'T_all'],as.integer(t[,paste0(c('Tma'),'_',sc)]))
                                     ,1,function(x) sum(x,na.rm=T))
      
      
      tab@data[t$X,'both_QT'] <- apply(cbind(tab@data[t$X,'both_QT'],as.integer(apply(cbind(
        apply(t[,paste0(c('Qmi','Qzf'),'_',sc)],1,all),
        t[,paste0(c('Tma'),'_',sc)]
      ),1,any)))
      ,1,function(x) sum(x,na.rm=T))
      
      
      for(v in vars) tab@data[t$X,v] <- apply(cbind(tab@data[t$X,v],as.integer(t[,paste0(v,'_',sc)]))
                                              ,1,function(x) sum(x,na.rm=T))
      
    }
    
    if(filter_lentic_out) out <- paste0(dir_model,clmod,'/SR_tab_',sc,'_lotic.rds')
    if(filter_lotic_out) out <- paste0(dir_model,clmod,'/SR_tab_',sc,'_lentic.rds')
    
    saveRDS(tab,out)
    
  },comboscen,mc.cores = length(comboscen),SIMPLIFY = F)
  
  
  cat('\n\n..done!\n')
  cat(paste0(rep('-',30)),'\n\n')
  
}
