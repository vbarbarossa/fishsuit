### FUNCTIONS

### load packages invisible
libinv <- function(packages){
  
  invisible(
    lapply(packages, function(x) suppressPackageStartupMessages( require(x, character.only = T, quietly = T)))
  )
  
}


### create directory if it does not exist
dir_ <- function(name_dir){
  if(!dir.exists(name_dir)) dir.create(name_dir,recursive = T)
  return(name_dir)
}

### crop raster based on polygon
crop_by_poly <-  function(in_ras_path,out_ras_path,crop_poly_path){
  
  system(
    paste0(
      'gdalwarp --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline ',
      crop_poly_path,' ',in_ras_path,' ',out_ras_path,' -overwrite'
    )
  )
  
  return(out_ras_path)
}

occurrence_tab <- function(gradient,spp.ranges.tab,f){
  
  SR = unlist(lapply(gradient,
                     function(g) sum(apply(as.matrix(spp.ranges.tab[,c('min','max')]), 1, findInterval, x=g) == 1))
  )
  
  PNOF = 1 - SR/max(SR)
  rSR = SR/max(SR)
  
  return(data.frame(Q = gradient,SR = SR, rSR = rSR, PNOF = PNOF, fac = as.factor(f)))
  
} ### end function

crop_ <- function(src,dst,crop_layer) {
  
  if(!file.exists(dst)){
    system(
      paste0(
        'ogr2ogr -overwrite -clipdst ',crop_layer,' ',dst,' ',src
      )
    )
    
  }
  
}

#function to extract annual averages from monthly data in ncdf files
#x (string)= path to the ncdf file
#varname (string)= the variable of the ncdf to be extracted
#band_start_month (integer)= band of the ncdf file corresponding to the first month (january usually) to be considered in the averaging
#seq_years (integer)= sequence of years to be extracted, must be continuous e.g. 1971,1972,1973 not 1971,1974,1976
extract_average <- function(x,varname,band_start_month,seq_years){
  library(raster)
  #set length timespan
  n_years = length(seq_years)
  #set sequence of "Januaries"
  nseq = seq(band_start_month,((n_years-1)*12+band_start_month),12)
  #calculate average
  for(i in 1:n_years){
    sum = raster(x, varname=varname, band = nseq[i])
    for(j in 1:11){
      sum = sum + raster(x, varname=varname, band = (nseq[i]+j))
    }
    assign(paste0('mean_',seq_years[i]),sum/12)
  }
  
  mean = get(paste0('mean_',seq_years[1]))
  for(i in 2:n_years){
    mean = mean + get(paste0('mean_',seq_years[i]))
  }
  mean = mean/n_years
  
  return(mean)
}

match_path <- function(text_file,matching_string){
  return(
    gsub(" ","",
         unlist(strsplit(text_file[grep(matching_string,text_file)],'=' ))[2]
    )
  )
}
