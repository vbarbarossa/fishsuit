#sbatch --array=1-5

g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source('config.R'); # always load as functions are loaded within this script

libinv(c('raster','foreach','dplyr','sf'))

clmod = climate_models[g]
scen = scenarios[scenarios != 'hist'] #<
years = floor(mean(timespan_scen)) #<


cat(paste0(rep('-',30)),'\n')
cat('thresholds used:\n')
cat(paste(vars,parse=' '))
cat('\n')
cat(paste(thresholds,parse=' '))
cat('\n',paste0(rep('-',30)),'\n\n')

# make one table per GCM
niche_list <- lapply(
  seq_along(vars),
  function(j){
    
    d <- readRDS(paste0('proc/',clmod,'/niches/',vars[j],'.rds')) %>% 
      dplyr::select(id_no = ID,paste0(thresholds[j],'%')) %>% as_tibble()
    colnames(d)[2] <- paste0(vars[j],'_',thresholds[j],'%')
    
    return(d)
  }
)

# make sure all ids correspond properly
niche <- niche_list[[1]] %>% as_tibble()
for(i in seq_along(vars)[-1]){
  niche <- inner_join(niche,niche_list[[i]])
}

niche <- niche %>% 
  mutate(id_no = as.character(id_no)) %>%
  filter(id_no %in% read.csv('proc/thresholds_average_filtered.csv')$id_no)

# filter based on overall species from threshold_selection.R
write.csv(niche,paste0('proc/',clmod,'/niches/niches_filtered.csv'))

cat(paste0(rep('-',30)),'\n\n')

#----------------------------------------------------------------------------------------------------

cat(paste0(rep('-',30)),'\n')
cat('modeling HABITAT SUITABILITY of each species..\n')

# make tab niches
nich <- niche %>% as.data.frame()
colnames(nich) <- c('id_no',paste0(vars,'_th'))

dir_out <- dir_(paste0('proc/',clmod,'/modelled_occurrence//'))

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
# when need to build the table, can simply create three columns for scen, year and warmt 
# and assign categories by splitting the string

# function that takes in input a species and calculates the suitability 
# of each cell based on the thresholds
build_table <- function(n){
  
  sp <- nich[n,'id_no']
  
  pts <- readRDS(paste0('proc/ssp/single_points/',sp,'.rds'))
  
  tab <- cbind(
    #includes info on which cells are
    cbind(pts %>% as.data.frame() %>% dplyr::select(-geometry),st_coordinates(pts))[,c('area','X','Y')], #includes info on which cells are
    
    foreach(v = 1:length(vars),.combine = 'cbind') %do% {
      
      # first check hist data
      r <- raster(paste0('proc/',clmod,'/pcrglobwb_processed/merged/',
                         vars[v],'_hist.tif'))
      
      val <- extract(r,pts)
      
      if(vars[v] == 'Tma'){
        return(val <= nich[n,paste0(vars[v],'_th')])
      }
      if(vars[v] == 'Tmi'){
        return(val >= nich[n,paste0(vars[v],'_th')])
      }
      if(vars[v] == 'Qma'){
        return(val <= nich[n,paste0(vars[v],'_th')])
      }
      if(vars[v] == 'Qzf'){
        return(val <= nich[n,paste0(vars[v],'_th')])
      }
      if(vars[v] == 'Qmi'){
        return(val >= nich[n,paste0(vars[v],'_th')])
      }
    },
    foreach(v = 1:length(vars),.combine = 'cbind') %do% {
      
      foreach(s = 1:length(comboscen),.combine = 'cbind') %do% {
        
        r <- raster(paste0('proc/',clmod,'/pcrglobwb_processed/merged/',
                           vars[v],'_',comboscen[s],'.tif'))
        
        val <- extract(r,pts)
        
        if(vars[v] == 'Tma'){
          return(val <= nich[n,paste0(vars[v],'_th')])
        }
        if(vars[v] == 'Tmi'){
          return(val >= nich[n,paste0(vars[v],'_th')])
        }
        if(vars[v] == 'Qma'){
          return(val <= nich[n,paste0(vars[v],'_th')])
        }
        if(vars[v] == 'Qzf'){
          return(val <= nich[n,paste0(vars[v],'_th')])
        }
        if(vars[v] == 'Qmi'){
          return(val >= nich[n,paste0(vars[v],'_th')])
        }
        
      }
      
    }
  )
  #assign col names to the table
  colnames(tab) <- c('area','x','y',
                     paste0(vars,'_hist'),
                     paste(rep(vars,each=length(comboscen)),comboscen,sep='_')
  )
  
  # maybe better write as RDS?
  # write.csv(tab,paste0(dir_out,sp,'.csv'),row.names = T) #impo to keep number of rows for later modelling
  saveRDS(tab,paste0(dir_out,sp,'.rds'))
  
}

invisible(
  parallel::mcmapply(build_table,1:nrow(nich),mc.silent = TRUE,mc.cores = ncores)
)

cat('..ended!\n')
cat(paste0(rep('-',30)),'\n\n')

#-------------------------------------------------------------------------------------------------------
# ESH_tab

cat(paste0(rep('-',30)),'\n')
cat('Build ESH_tab\n')

tabs <- list.files(dir_out,pattern='.rds')
ids <- as.character(
  sapply(
    list.files(dir_out,pattern='.rds'),function(x) strsplit(strsplit(x,'_')[[1]][1],'\\.')[[1]][1]
  )
)


# function to calculate ESH for each species-tab
esh_calc <- function(sp){ #,sn
  # read the table
  t <- readRDS(paste0(dir_out,sp,'.rds'))
  # <<<<<<<<<<<<< should also load the other tabl to et cell area, or select from template (but it's bigger)
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

out <- paste0(dir_model,clmod,'/ESH_tab.csv')
write.csv(ESH_tab,out,row.names = F)

cat('successfully written to ',out,'\n')
cat(paste0(rep('-',30)),'\n\n')


#----------------------------------------------------------------------------------------------------
# SUMMARIZE at the cell-level

cat(paste0(rep('-',30)),'\n')
cat('Build SR_tab\n')

for(sc in comboscen){
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
  
  saveRDS(tab,paste0(dir_model,clmod,'/SR_tab_',sc,'.rds'))
  
}


cat('\n\n..done!\n')
cat(paste0(rep('-',30)),'\n\n')



