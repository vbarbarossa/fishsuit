# RU CLUSTER --> could not install rredlist on cartesius

# retrieve habitat type
library("rredlist"); library(foreach);
source('R/functions.R'); source('R/_dir_loc.R');

dir_habitat <- dir_(paste0(dir_proc,'ssp/habitats/'))
dir_iucn_data <- paste0(dir_data,'IUCN/FW_FISH_20181113/')

#my personal token
token <- 'd361026f05b472e57b0ffe1fa5c9a768aaf3d8391abbb464293e9efe2bbbf733'

#iucn fish data
dat <- rbind(
  foreign::read.dbf(paste0(dir_iucn_data,'FW_FISH_PART_1.dbf')),
  foreign::read.dbf(paste0(dir_iucn_data,'FW_FISH_PART_2.dbf'))
)

# get unique species names
sp_names <- unique(as.character(dat[,'binomial']))

#in case the scripts fails, do only remaining species
sp_names_completed <- as.character(sapply(list.files(dir_habitat),function(x) strsplit(x,'\\.')[[1]][1]))
sp_names <- sp_names[!(sp_names %in% sp_names_completed)] 

# retrieve the habitat type per species
for(i in sp_names) {
  
  hab <- rl_habitats(as.character(i),key = token)$result$habitat
  saveRDS(hab,paste0(dir_habitat,i,'.rds'))
  
  Sys.sleep(2)
  
}
