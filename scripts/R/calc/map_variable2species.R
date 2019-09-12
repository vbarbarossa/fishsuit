#Valerio Barbarossa, 12 Sep 2019
# function that takes in input:
# -an environmental variable (raster layer)
# -species occurrence range (points) 
# and outputs:
# -table with the distribution of the values within the range (quantiles)

map_variable2species <- function(
  # define function variables
  infile_lyr #.tif
  ,infiles_pts #.rds # path to the point shapefiles
  ,ids # vector of ID names
  ,outfile_name # name of output table (NO EXTENTION). The file is saved in csv and rds formats
  ,NC = 1 # number of cores to use
  ,dir_tmp = dir_(paste0('tmp',sample(1:10**5,1))) # set temporary directory for calculations, if not specified a random one in used
){
  
  #load packages
  library('valerioUtils')
  libinv(c('raster','foreach','parallel'))
  
  # set raster package tmp dir
  rasterOptions(tmpdir = dir_tmp)
  
  # load variable raster layer in memory
  lyr <- raster(infile_lyr)
  
  # extract lyr values to infiles_pts in parallel -----------------------------------------------
  t_range <- do.call(
    'rbind',
    mcmapply(function(x){
      # load points in memory
      pts <- readRDS(x)
      
      # extract values
      val <- extract(lyr,pts)
      # exclude eventual NAs
      val <- val[!is.na(val)]
      
      # retrieve quantiles of the values distribution
      q <- quantile(val,seq(0,1,0.005))
      
      # return the data as data.frame
      return(
        cbind(data.frame(no.grids = length(val), mean = mean(val), sd = sd(val)),
              t(as.data.frame(q)))
      )
      
    },infiles_pts,SIMPLIFY = F,mc.cores = NC)
  )
  row.names(t_range) <- NULL
  t_range <- cbind(data.frame(ID = ids),t_range)
  #---------------------------------------------------------------------------------------------
  
  # save results
  write.csv(t_range,paste0(outfile_name,'.csv'),row.names = FALSE)
  saveRDS(t_range,paste0(outfile_name,'.rds'))
  
  # remove temporary directory
  rm_(dir_tmp)
  
}
