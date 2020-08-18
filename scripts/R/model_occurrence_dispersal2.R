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

dir_out <- dir_(paste0('proc/',clmod,'/modelled_occurrence_dispersal2/'))

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
  
  pts <- readRDS(paste0('proc/ssp/single_points_dispersal/',sp,'.rds')) %>% dplyr::select(-ws)
  
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
  
  #<<<<<<<<<<< version2 modification
  # add a column that says whether the row/cell belongs to the original occurrence range (1)
  # or is part of the expanded part (0)
  tab$occ <- 0
  tab$occ[row.names(tab) %in% as.character(readRDS(paste0('proc/ssp/single_points/',sp,'.rds'))$row_no)] <- 1
  #<<<<<<<<<<<
  
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


# function to calculate RT for each species-scenario combination
esh_calc <- function(sp){ #,sn
  
  # function that calculater % of Range Threatened (RT) as
  # RT = L/(A+Ae-Le)
  # where 
  # L is the area of range lost within the original range
  # A is the area of the original range
  # Ae is the area of the expanded part of the exapnded range
  # Le is the area of range lost within the expanded part of the expanded range
  calc_RT <- function(va){
    if(length(va) > 1){
      return(
        sum(t[(!apply(t[,va],1,all)) & (t$occ == 1),'area'])/(d$area_adjusted_km2[y] - sum(t[(!apply(t[,va],1,all)) & (t$occ == 0),'area'])) *100
      )
    }else{
      return(
        sum(t[(!t[,va]) & (t$occ == 1),'area'])/(d$area_adjusted_km2[y] - sum(t[(!t[,va]) & (t$occ == 0),'area'])) *100
      )
    }
    
  }
  
  # read the table
  t <- readRDS(paste0(dir_out,sp,'.rds'))
  t <- t[complete.cases(t),] #remove rows that contain one or more NAs
  
  d <- data.frame(
    id_no = rep(sp,length(comboscen)),
    no_cells = nrow(t),
    area_total_km2 = sum(t$area),
    ESH_all = NA,
    ESH_Q = NA,
    ESH_T = NA,
    ESH_QnT = NA,
    comboscen = comboscen #sn is the looped scen
  )
  for(v in vars) d[,paste0('ESH_',v)] <- NA
  
  for(y in 1:nrow(d)){
    
    sn = d$comboscen[y]
    
    all_vars_hist <- paste0(vars,'_hist')
    all_vars_scen <- paste0(vars,'_',sn)
    
    #filter out rows of ts that are already excluded in th due to the quantile threshold
    t <- t[which(apply(t[,all_vars_hist],1,all)),]
    
    # compute adjusted area after corrections
    d$area_adjusted_km2[y] <- sum(t$area)
    
    # need to count the number of false and calculate the associated area per grid-cell
    d$ESH_all[y] <- calc_RT(all_vars_scen)
    d$ESH_Q[y] <- calc_RT(paste0(c('Qmi','Qzf','Qma'),'_',sn))
    d$ESH_T[y] <- calc_RT(paste0(c('Tmi','Tma'),'_',sn))
    
    d$ESH_QnT[y] <- sum(t[!apply(
      cbind(
        apply(t[,paste0(c('Qmi','Qzf','Qma'),'_',sn)],1,all),
        apply(t[,paste0(c('Tmi','Tma'),'_',sn)],1,all))
      ,1,any) & (t$occ == 1),'area'])/(d$area_adjusted_km2[y] - sum(t[!apply(
        cbind(
          apply(t[,paste0(c('Qmi','Qzf','Qma'),'_',sn)],1,all),
          apply(t[,paste0(c('Tmi','Tma'),'_',sn)],1,all))
        ,1,any) & (t$occ == 0),'area']))*100
    
    for(v in vars) d[y,paste0('ESH_',v)] <- calc_RT(paste0(v,'_',sn))
    
    
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

out <- paste0('proc/',clmod,'/ESH_tab_dispersal2.csv')
write.csv(ESH_tab,out,row.names = F)

cat('successfully written to ',out,'\n')
cat(paste0(rep('-',30)),'\n\n')


#----------------------------------------------------------------------------------------------------
# SUMMARIZE at the cell-level

cat(paste0(rep('-',30)),'\n')
cat('Build SR_tab\n')

# in the dispersal scenario for each variable assessed 
# also the number of species for the historical scenario need to be calculated
# to be able to exclude in the present expanded ranges
# those cells that will be unsuitable in the future (within the expanded part)
# the below function calculates that
calc_occ <- function(d,va,occ_m=FALSE,any_m = FALSE, both_m = FALSE){
  
  if(length(va) > 1){
    d$o <- as.integer(apply(t[,paste0(va,'_',sc)],1,all))
    if(any_m) d$o <- as.integer(apply(t[,paste0(va,'_',sc)],1,any))
    if(both_m){
      d$o <- as.integer(apply(cbind(
        apply(t[,paste0(c('Qmi','Qzf','Qma'),'_',sc)],1,all),
        apply(t[,paste0(c('Tmi','Tma'),'_',sc)],1,all)
      ),1,any))
    }
    if(occ_m) d$o[d$occ == 1] <- 1
    return(d$o)
  }else{
    d$o <- as.integer(t[,paste0(va,'_',sc)])
    if(occ_m) d$o[d$occ == 1] <- 1
    return(d$o)
  }
  
}

for(sc in comboscen){
  
  tab <- readRDS('proc/ssp/points_template.rds') %>% as.data.frame() %>% dplyr::select(-geometry)
  for(v in c(vars,'Q_all','T_all','any','all','both_QT','occ')) tab[,v] <- NA
  for(v in colnames(tab)[-c(1,2,ncol(tab))]) tab[,paste0('occ_',v)] <- NA
  
  
  for(i in ids){
    
    t <- readRDS(paste0(dir_out,i,'.rds'))
    t$X <- as.integer(row.names(t)) #the way it was coded based on csv
    th <- t[,paste0(vars,'_hist')] #historical
    
    t <- t[which(apply(th,1,all)),]
    
    # this will not be used for the final calculatiopns, it is simply the richness with fully expanded ranges
    tab$occ[t$X] <- apply(cbind(tab$occ[t$X],rep(1,nrow(t))),1,function(x) sum(x,na.rm=T)) 
    
    # all
    tab[t$X,'all'] <- apply(
      cbind(
        tab[t$X,'all'],
        calc_occ(t,c('Qmi','Qzf','Qma','Tma','Tmi'))
      )
      ,1,function(x) sum(x,na.rm=T))
    
    tab[t$X,'occ_all'] <- apply(
      cbind(
        tab[t$X,'occ_all'],
        calc_occ(t,c('Qmi','Qzf','Qma','Tma','Tmi'),occ_m = T)
      )
      ,1,function(x) sum(x,na.rm=T))
    
    
    tab[t$X,'any'] <- apply(
      cbind(
        tab[t$X,'any'],
        calc_occ(t,c('Qmi','Qzf','Qma','Tma','Tmi'),any_m = T)
      )
      ,1,function(x) sum(x,na.rm=T))
    
    tab[t$X,'occ_any'] <- apply(
      cbind(
        tab[t$X,'occ_any'],
        calc_occ(t,c('Qmi','Qzf','Qma','Tma','Tmi'),occ_m = T, any_m = T)
      )
      ,1,function(x) sum(x,na.rm=T))
    
    
    tab[t$X,'Q_all'] <- apply(
      cbind(
        tab[t$X,'Q_all'],
        calc_occ(t,c('Qmi','Qzf','Qma'))
      )
      ,1,function(x) sum(x,na.rm=T))
    
    tab[t$X,'occ_Q_all'] <- apply(
      cbind(
        tab[t$X,'occ_Q_all'],
        calc_occ(t,c('Qmi','Qzf','Qma'),occ_m = T)
      )
      ,1,function(x) sum(x,na.rm=T))
    
    
    tab[t$X,'T_all'] <- apply(
      cbind(
        tab[t$X,'T_all'],
        calc_occ(t,c('Tmi','Tma'))
      )
      ,1,function(x) sum(x,na.rm=T))
    
    tab[t$X,'occ_T_all'] <- apply(
      cbind(
        tab[t$X,'occ_T_all'],
        calc_occ(t,c('Tmi','Tma'),occ_m = T)
      )
      ,1,function(x) sum(x,na.rm=T))
    
    tab[t$X,'both_QT'] <- apply(
      cbind(
        tab[t$X,'both_QT'],
        calc_occ(t,c('Qmi','Qzf','Qma','Tma','Tmi'),both_m = T)
      ),1,function(x) sum(x,na.rm=T))
    
    tab[t$X,'occ_both_QT'] <- apply(
      cbind(
        tab[t$X,'occ_both_QT'],
        calc_occ(t,c('Qmi','Qzf','Qma','Tma','Tmi'),both_m = T,occ_m = T)
      ),1,function(x) sum(x,na.rm=T))
    
    for(v in vars) tab[t$X,v] <- apply(
      cbind(
        tab[t$X,v],
        calc_occ(t,v)
      ),1,function(x) sum(x,na.rm=T))
    
    for(v in vars) tab[t$X,paste0('occ_',v)] <- apply(
      cbind(
        tab[t$X,paste0('occ_',v)],
        calc_occ(t,v,occ_m = T)
      ),1,function(x) sum(x,na.rm=T))
    
    
  }
  
  saveRDS(tab,paste0('proc/',clmod,'/SR_tab_dispersal2_',sc,'.rds'))
  
}


cat('\n\n..done!\n')
cat(paste0(rep('-',30)),'\n\n')

