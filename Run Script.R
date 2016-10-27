#Install required packages
packages <- c("PBSmapping", "sp", "rgeos", "maptools", "raster","SDMTools",
              "car","nlme", "multcomp","lsmeans","MuMIn","lubridate","adehabitatHR",
              "geometry", "fpc", "FNN","circular")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}


