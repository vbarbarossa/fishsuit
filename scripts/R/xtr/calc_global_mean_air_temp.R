library(raster)

dir_input <- '/lustre5/0/milkun/Input_for_water_temperature/'
dir_out <- 'global_Ta/'
dir.create(dir_out)

climate_models <- c('GFDL','HadGEM','IPSL','MIROC','NorESM')
complement <- c('-ESM2M','2-ES','-CM5A-LR','-ESM-CHEM','1-M')

scenarios <- c('2p6','4p5','6p0','8p5')

for(i in 1:length(climate_models)){
  
  # historical 1951:2005
  yrs <- 1953:2003
  
  vals <- data.frame(year = yrs, Tmean = NA)
  
  for(y in 1:length(yrs)){
    
    r <- raster(paste0(dir_input,climate_models[i],'/',climate_models[i],complement[i],
                       '_annualT_hist.nc'),band=y)
    
    vals[y,'Tmean'] <- mean(values(r),na.rm=T) - 273.15
    
  }
  
  write.csv(vals,paste0(dir_out,climate_models[i],'_hist.csv'),row.names = F)
  
  
  
  
  # future 2006:2099
  yrs <- 2008:2097
  
  for(scen in scenarios){
    valsf <- data.frame(year = yrs, Tmean = NA)
  
    for(y in 1:length(yrs)){
      
      r <- raster(paste0(dir_input,climate_models[i],'/',climate_models[i],complement[i],
                         '_annualT_',scen,'.nc'),band=y)
      
      valsf[i,'Tmean'] <- mean(values(r),na.rm=T) - 273.15
      
    }
    
    write.csv(valsf,paste0(dir_out,climate_models[i],'_',scen,'.csv'),row.names = F)
    
  }
  
  
}

