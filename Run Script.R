#Install required packages
packages <- c("PBSmapping", "sp", "rgeos", "maptools", "raster","SDMTools",
              "car","nlme", "multcomp","lsmeans","MuMIn","lubridate","adehabitatHR",
              "geometry", "fpc", "FNN","circular","rgdal","ggmap","ggplot2","mapplots",
              "lmerTest","zoo","AICcmodavg")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(devtools)
install_github("yihui/printr")

