##Clear R Workspace
rm(list=ls())

#Required R Packages:
require(PBSmapping)
require(sp)
require(rgeos)
require(maptools)
require(raster)

#load functions
source('Functions/VPS Error Functions.R')

#####Step 1#### 
#Inputting and Importing Data

#VPS data 2013
#Load VPS csv
all_data1<-read.csv("Data/VPS 2012/ALL-CALC-POSITIONS.csv")
all_data2<-read.csv("Data/VPS 2013/ALL-CALC-POSITIONS1.csv")
all_data3<-read.csv("Data/VPS 2013/ALL-CALC-POSITIONS2.csv")
all_data4<-read.csv("Data/VPS 2013/ALL-CALC-POSITIONS3.csv")
all_data5<-read.csv("Data/VPS 2014/ALL-CALC-POSITIONS.csv")

#combine summer and winter datasets
all_data<-rbind(all_data,all_data2,all_data3)
#remove winter data set
rm(all_data2,all_data3)



#load shoreline polygon and projection
shore_outline<-readShapeSpatial("~/Documents/Masters Thesis/Data/GIS/qgis/Alexie.shp")
proj4string(shore_outline)<-CRS("+proj=utm +zone=11 +ellps=WGS84")



#Prep data
#assign timezone for VPS data and project in UTM
data<-load_coords(data=all_data, output_timezone="MST")

rm(all_data)

#extract data during lake trout spawn 2013
data<-data[data$DATETIME>=as.POSIXct("2013-05-15",tz="MST") & 
           data$DATETIME<as.POSIXct("2013-12-15",tz="MST"),]

#FOCAL FISH POSITIONS
tag<-c("LT-31","LT-32","LT-33","LT-34","LT-35","LT-37","LT-38","LT-40","LT-41","LT-43","LT-44")

pos <- data[grep(paste(tag,collapse="|"), data$TRANSMITTER), ]
pos$TRANSMITTER<-droplevels(pos$TRANSMITTER)


rm(data)
