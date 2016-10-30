##Clear R Workspace
#rm(list=ls())

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


#Load VPS csv
#VPS data 2012
all_data1<-read.csv("~/Data/VPS 2012/ALL-CALC-POSITIONS.csv")
all_data2<-read.csv("~/Data/VPS 2012/ALL-CALC-POSITIONS1.csv")
all_data3<-read.csv("~/Data/VPS 2012/ALL-CALC-POSITIONS2.csv")
all_data4<-read.csv("~/Data/VPS 2012/ALL-CALC-POSITIONS3.csv")
#VPS data 2013
all_data5<-read.csv("~/Data/VPS 2013/ALL-CALC-POSITIONS1.csv")
all_data6<-read.csv("~/Data/VPS 2013/ALL-CALC-POSITIONS2.csv")
all_data7<-read.csv("~/Data/VPS 2013/ALL-CALC-POSITIONS3.csv")
#VPS data 2014
all_data8<-read.csv("~/Data/VPS 2014/ALL-CALC-POSITIONS.csv")

#combine summer and winter datasets
all_data<-rbind(all_data1,all_data2[,c(2:17)],all_data3,all_data4,all_data5,
                all_data6,all_data7,all_data8)
#remove winter data set
rm(all_data1,all_data2,all_data3,all_data4,all_data5,all_data6,all_data7,all_data8)



#load shoreline polygon and projection
shore_outline<-readShapeSpatial("~//Data/GIS/Alexie.shp")
proj4string(shore_outline)<-CRS("+proj=utm +zone=11 +ellps=WGS84")



#Prep data
#assign timezone for VPS data and project in UTM
data<-load_coords(data=all_data, output_timezone="MST")

rm(all_data)
#Sort detetctions by datetime
data<-data[order(data$DATETIME),]

#extract data during lake trout spawn 2012-2014 spawn season
data<-data[data$DATETIME>=as.POSIXct("2012-08-15",tz="MST") & 
           data$DATETIME<as.POSIXct("2012-12-15",tz="MST") |
           data$DATETIME>=as.POSIXct("2013-08-15",tz="MST") & 
           data$DATETIME<as.POSIXct("2013-12-15",tz="MST") |
           data$DATETIME>=as.POSIXct("2014-08-15",tz="MST") & 
           data$DATETIME<as.POSIXct("2014-12-15",tz="MST")
           ,]

#extract lake trout positions only
pos <- data[grep("LT", data$TRANSMITTER), ]

pos$TRANSMITTER<-droplevels(pos$TRANSMITTER)


rm(data)
