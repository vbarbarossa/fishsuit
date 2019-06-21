library(raster)
source('R/_dir_loc.R'); source('R/functions.R')

dir_points <- dir_(paste0(dir_proc,'ssp/ssp_points/'))

# read template file and ids
ids <- sapply(strsplit(list.files(dir_points,pattern = '.rds'),'\\.'),function(x) x[1])
tab <- readRDS(paste0(dir_proc,'ssp/tab_template.rds'))

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# here needs to read the actual ids used after filtering for the study

index=1
pb <- txtProgressBar(min = 1, max = length(ids), style = 3)

for(i in ids){
  
  t <- readRDS(paste0(dir_points,i,'.rds'))
  t@data$X <- as.numeric(row.names(t@data))
  t <- t@data
  #all
  tab@data$occ[t$X] <- apply(cbind(tab@data$occ[t$X],rep(1,nrow(t))),1,function(x) sum(x,na.rm=T)) 
  
  setTxtProgressBar(pb, index)
  index=index+1
  
}

close(pb)

tab_filtered <- tab[!is.na(tab@data$occ),]

ras <- rasterize(tab_filtered,raster(resolution=1/12))

saveRDS(ras,paste0(dir_proc,'SR_raster_original.rds'))
