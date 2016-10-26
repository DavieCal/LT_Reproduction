###############################
#DATA ERROR SCRIPT
###############################

#Author: Dave Callaghan
#Date Created:07-04-2014
#Date Edited:05-12-2014
#
###############################
#Libraries
library(chron)
library(rgdal)
library(vgma)
library(rgeos)
library(maptools)
###############################
#Set Variables



setwd("~/Documents/Masters Thesis/Data/VPS")
AllData<-read.csv("VPS_UTM_ATTR_2012-2013.csv")
AllData<-data

#Set working directory
setwd("~/Documents/R-Directory/AniMove")

#Load all data and project
source("LoadData.R")

#Save current parsettings
opar <- par() 

#create subsets

start=chron(dates = "06/12/13", times = "0:00:00",
             format = c(dates = "m/d/y", times= "h:m:s"))

end=chron(dates = "06/30/13", times = "0:00:00",
           format = c(dates = "m/d/y", times= "h:m:s"))
#subset by data
sub<-data[data$timestamp>start & data$timestamp<end,]

#determine available transmitters for subset
unique(sub$TRANSMITTER)

#32, 39, 34, 44, 40, 37, 43, 42, 41
tag=41
f=9
#subset by tag
data.unclean<-sub[sub$TRANSMITTER==paste("LT-",tag,sep=""),]



setwd("~/Documents/R-Directory/AniMove")
#remove duplicates from data
###NOTE: Need to add options for data cleaning (ie which duplicate to keep)
source("CleanData.R")

points(catch_data$UTM.X[f],catch_data$UTM.Y[f],col="red")
#sample
assign(paste("tag",tag,sep=""),data.clean)

#32, 39, 34, 44, 40, 37, 43, 42, 41
tag=41

#setwd("~/Documents/R-Directory/Masters Thesis/DATA/Alexie/Fish Data/")
setwd("~/Documents/R-Directory/homing")
save(list=paste("tag",tag,sep=""),
     file=paste("image_tag",tag,".rdata",sep=""))
#rm(list=paste("tag",tag,sep=""))

fem<-c("image_tag32.rdata","image_tag35.rdata","image_tag38.rdata",
            "image_tag41.rdata","image_tag43.rdata","image_tag44.rdata")
load("image_tag32.rdata")
load("image_tag35.rdata")
load("image_tag38.rdata")
load("image_tag41.rdata")
load("image_tag43.rdata")
load("image_tag44.rdata")

#output assignment
samp<-data.clean

###############################
# Start processing
print("Set variables and start processing")
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")
 

depth.vario<-getVariogram(samp$timestamp,samp$DEPTH,samp$DEPTH)
#all data
plot(samp.vario,tc=1,ylim=c(0,175000),xlim=c(0,20))
plot(samp.vario,tc = 1,xlim=c(0,0.1),ylim=c(0,40))

#best data
#samp2.vario<-getVariogram(samp$timestamp,samp$UTM.X,samp$UTM.Y)
plot(samp2.vario,tc=1,ylim=c(0,175000),xlim=c(0,20))
plot(samp2.vario,tc = 1,xlim=c(0,0.01),ylim=c(0,40))

print("Set variables and start processing")
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

#good data
#samp3.vario<-getVariogram(samp$timestamp,samp$UTM.X,samp$UTM.Y)
plot(samp3.vario,tc=1,ylim=c(0,175000),xlim=c(0,20))
plot(samp3.vario,tc = 1,xlim=c(0,0.01),ylim=c(0,40))

# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")

save(dbbmm_stack ,file="dbbmm_stack.rdata")

setwd("~/Documents/R-Directory/AniMove")
load("image_DataError_april08.rdata")

tag1<-data.clean
tag3<-data.clean
tag22<-data.clean
tag31<-data.clean
