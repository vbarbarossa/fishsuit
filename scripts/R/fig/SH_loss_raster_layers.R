source('R/4targets/MASTER.R');

library(raster); library(foreach); library(sf)
#> ENSEMBLE RASTER LAYERS -------------------------------------------------------------------------------------

# attribute table with year/warming target/scenario/climate model combos
library(foreach)
attr.tab <- foreach(i = seq_along(warming_targets),.combine = 'rbind') %do% {
  data.frame(
    clmod = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[1]),
    scen = sapply(strsplit(as.character(warming_tab[,1]),'_'),function(x) x[2]),
    year = warming_tab[,i+1],
    warmt = warming_targets[i])
}
attr.tab <- attr.tab[!is.na(attr.tab$year),]
row.names(attr.tab) <- NULL

# function that takes in input the row of the attr.tab and gives the rasterized relative losses per grid cell
rasterize_rel_losses <- function(x,dir_root = 'fishsuit_completeRun_warming_4targets'){
  tab <- readRDS(paste0(dir_root,'/',atab$clmod[x],'/SR_tab_',
                        atab$scen[x],'_',atab$warmt[x],'C_',atab$year[x],'.rds'))
  tab <- tab[!is.na(tab@data$occ),]
  tab@data <- data.frame(rel_all = (tab@data$occ - tab@data$all)/tab@data$occ)
  return(
    rasterize(tab,raster(ext=extent(tab),resolution=1/12),field = 'rel_all')
  )
}

list.rasters.name <- paste0(dir_mod,dir_('figs/data/'),'raster_layers_figure_1.rds')

if(file.exists(list.rasters.name)){
  list.rasters <- readRDS(list.rasters.name)
}else{
  list.rasters <- foreach(n = seq_along(warming_targets)) %do% {
    
    atab <- attr.tab[attr.tab$warmt == warming_targets[n],]
    
    ncores = 7
    # parallelized ---
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    br <- brick(
      foreach(i = 1:nrow(atab),.packages='raster') %dopar% rasterize_rel_losses(x=i)
    )
    parallel::stopCluster(cl)
    # ---
    return(
      calc(br,function(x) median(x,na.rm=T))
    )
    
  }
  names(list.rasters) <- warming_targets
  
  saveRDS(file = list.rasters.name,object = list.rasters)
  
}

#> MAPS -------------------------------------------------------------------------------------------------------

for(i in 1:length(list.rasters)){
  
  r <- list.rasters[[i]]
  rbase <- raster('ldd.asc')
  rbase[!is.na(rbase)] <- -1
  
  res(rbase) <- res(r)
  rbase <- crop(rbase,extent(r))
  extent(rbase) <- extent(r)
  
  r[is.na(r) & rbase==-1] <- NA # cells with no species data available
  
  writeRaster(r,paste0(dir_('fishsuit_completeRun_warming_4targets/figshare/'),'cumulative_range_loss_',names(list.rasters)[i],'.tif'),
              format = 'GTiff',overwrite = T)
  
  
}

