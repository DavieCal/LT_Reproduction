##################
#Data Cleaning
#
#removes duplicate time stamps and positions.
#

#NOTE: setup rules for selecting same timestamps...
#      alternatively do manually.
#data_unclean= VPS data
#min_interval= minimum ping interval
#velocity_cut= max velocity of animal

#####read in bathymetric raster#####    
#if the raster already exists it will skip this step
if(!exists("r")){
  require(SDMTools)
  bathy<-read.asc("~/Data/GIS//bathy.asc")
  #assign raster class object
  r<-raster.from.asc(bathy)
  #assign projection
  projection(r)<- "+proj=utm +zone=11 +ellps=WGS84"
  r<-r*-1
  # unload, some functions conflict with raster functions
  detach("package:SDMTools", unload=TRUE)
  
}

require(rgeos)
require(maptools)
if(!exists("shore_outline")){
  #load shoreline polygon and projection
  shore_outline<-readShapeSpatial("~/Data/GIS/Alexie.shp")
  proj4string(shore_outline)<-CRS("+proj=utm +zone=11 +ellps=WGS84")
}
##test data
#data<-data_unclean[[1]]
#min_interval=80
#velocity_cut=1.2

bioFilter<-function(data=data,min_interval=min_interval,velocity_cut=velocity_cut){

#order data by data
data<-data[order(data$DATETIME), ]


#reports any duplicates for a single individual in terms of timestamps
if (any(duplicated(unclass(data$DATETIME)))){
  print(paste("duplicate timestamps:", sum(duplicated(unclass(data$DATETIME)))))
  #readline(prompt="Press [enter] to continue")
  data.clean <- data[!duplicated(data[,c("DATETIME", "UTM.X", "UTM.Y")]),]
}else {
  print("No duplicates")
  data.clean<-data
}

#####remove points less than minimum interval

#keep all positions >= than minimum interval seconds apart and the first point (NA time diff)


data.clean<-data.clean[order(data.clean$DATETIME),]
data.clean$TIMEDIFF<-NA
data.clean$TIMEDIFF[2:dim(data.clean)[1]]<-with(data.clean,difftime(DATETIME[-1],DATETIME[-length(DATETIME)],unit="secs")) 
data.clean1<-data.clean[is.na(data.clean$TIMEDIFF) | data.clean$TIMEDIFF>=min_interval,] 
#continue removing until no positions < minimum time apart
while(!dim(data.clean1)[1]== dim(data.clean)[1]){
  data.clean=data.clean1
  data.clean$TIMEDIFF<-NA
  data.clean$TIMEDIFF[2:dim(data.clean)[1]]<-with(data.clean,difftime(DATETIME[-1],DATETIME[-length(DATETIME)],unit="secs")) 
  data.clean1<-data.clean[is.na(data.clean$TIMEDIFF) | data.clean$TIMEDIFF>=min_interval,]
}                                                                                                                                                         
rm(data.clean1)
cat("Data points removed: ",dim(data)[1]-dim(data.clean)[1],"\n")
##################
    
# Now its time to create the classes
# Classifying points based on criteria
# Criteria:
# 1) Point within shoreline or within 3 m of shoreline = 1
# 2) Position above bottom of lake or within 2 m = 1
# 3) Realistic swimming speed = 1
# 

##################
#In shoreline or with 3 m

require(rgeos)
require(maptools)
if(!exists("shore_outline")){
#load shoreline polygon and projection
shore_outline<-readShapeSpatial("~/Data/GIS/Alexie.shp")
proj4string(shore_outline)<-CRS("+proj=utm +zone=11 +ellps=WGS84")
}

#shore distance

if(!exists("shore.dist")){
  require(SDMTools)
  shore.dist<-read.asc("~/Data/GIS/shoredist2.asc")
  #assign raster class object
  shore.dist<-raster.from.asc(shore.dist)
  #assign projection
  projection(shore.dist)<- "+proj=utm +zone=11 +ellps=WGS84"
  detach("package:SDMTools", unload=TRUE)
}


#create spatial data frame and assign projection
coordinates(data.clean)<-~UTM.X+UTM.Y
proj4string(data.clean)<-CRS("+proj=utm +zone=11 +ellps=WGS84")

#extract distance from shore for each point
data.clean$SHOREDIST<-extract(shore.dist, data.clean)

inLake<-gContains(shore_outline,data.clean,byid=TRUE)
data.clean$SHOREDIST2<-ifelse(inLake[,1],data.clean$SHOREDIST*-1,data.clean$SHOREDIST)

rm(inLake)
#determine if points in polygon or with 3 m of the edge assign a value of 1 if TRUE
data.clean$INLAKE[data.clean$SHOREDIST2<=3]=1
#remove NA's
data.clean$INLAKE[is.na(data.clean$INLAKE)]=0


##################
#depth above bottom of lake


#####read in bathymetric raster#####    
#if the raster already exists it will skip this step
if(!exists("r")){
  require(SDMTools)
  bathy<-read.asc("~/Data/GIS/bathy.asc")
  #assign raster class object
  r<-raster.from.asc(bathy)
  #assign projection
  projection(r)<- "+proj=utm +zone=11 +ellps=WGS84"
  r<-r*-1
  # unload, some functions conflict with raster functions
  detach("package:SDMTools", unload=TRUE)
  
}


#extract depth at position and depth of bathymetry at position
bottom <- data.frame(data.clean$DEPTH,
                         extract(r, data.clean))
#determine difference between fish depth and possible depth
bottom$Diff<-bottom[,2]-bottom[,1]
#create above sedement column
bottom$AboveSed=0
#any difference of greater than -2 meters is give a 1
bottom$AboveSed[bottom$Diff>=-2]=1


data.clean$ABOVEBOTTOM=bottom$AboveSed


##################
#Realistic swim speed

#unclass coordinates and create a dataframe
data.unclass<-as.data.frame(coordinates(data.clean))
x<-data.unclass$UTM.X
y<-data.unclass$UTM.Y
#add x and y data as real and imaginary numbers in the z vector
z<-x+1i*y  

#determine step lengths
dz<-diff(z)
#steplength
s<-Mod(dz)
summary(s)

par(mfrow=c(1,3))
#Distribution of step lengths
hist(s, col="grey", bor="darkgrey", freq=FALSE,
     main="Histogram of steps lengths",
     xlab="meters")
lines(density(s), col=2, lwd=2)


t2<-data.clean$TIMEDIFF[-1]

hist(t2, col="grey", bor="darkgrey", freq=FALSE,
     main="Histogram of time steps",
     xlab="seconds")
lines(density(t2), col=2, lwd=2)

print("Time step summary")
summary(t2)


#calculate speed
v<-s/t2


hist(v, col="grey", bor="darkgrey", freq=FALSE,
     main="Histogram of speed steps",
     xlab="m/s")
lines(density(v), col=2, lwd=2)

print("Speed step summary")
summary(v)

#create speed column
data.clean$GROUNDSPEED=0
#assign values to destination
data.clean$GROUNDSPEED[2:length(data.clean)]=v

#create realistic swim column
data.clean$REALSPEED=0
#assign a value of 1 to swim speeds of less than or equal to velocity cut in m/s
data.clean$REALSPEED[data.clean$GROUNDSPEED<=velocity_cut]=1
#remove false positive unrealistic speeds(keep point after error point)
for (i in 2:dim(data.clean)[1]){
  if(data.clean$REALSPEED[i]==0 & data.clean$REALSPEED[i-1]==0){
    data.unclass<-as.data.frame(coordinates(data.clean[c(i-2,i),]))
    x<-data.unclass$UTM.X
    y<-data.unclass$UTM.Y
    #add x and y data as real and imaginary numbers in the z vector
    z<-x+1i*y  
    
    #determine step lengths
    dz<-diff(z)
    #steplength
    s<-Mod(dz)
    t<-as.numeric(difftime(data.clean$DATETIME[i],data.clean$DATETIME[i-2],units="secs"))
    data.clean$GROUNDSPEED[i]<-s/t
  }
}
data.clean$REALSPEED[data.clean$GROUNDSPEED<=velocity_cut]=1


##################
#Data Quality Score
data.clean$DATASCORE <- (data.clean$INLAKE + data.clean$REALSPEED + data.clean$ABOVEBOTTOM)

hist(data.clean$DATASCORE, 
     main="Histogram of Data Quality", 
     xlab="Data Quality Score")

cat("tag ", paste(data.clean$TRANSMITTER[1]),"\n")
cat("Total Positions = ", dim(data.clean)[1],"\n")

cat("Data Score","\n",names(table(data.clean$DATASCORE)),"\n",table(data.clean$DATASCORE),"\n")

par(mfrow=c(1,1))

return(as.data.frame(data.clean))
}

