
#plots 6.15 x 4.76
source('Functions/st-dbscan function.R')
library(lubridate)
library(adehabitatHR)
#library(rgl)
#library(plot3D)
#library(plot3Drgl)
library(raster)
library(geometry)
library(fpc)
library(FNN)
library(rgeos)
library(circular)

#Classify males, females and unknowns
male<-c("LT-1","LT-3","LT-8","LT-22","LT-31","LT-33",
        "LT-34","LT-36","LT-37","LT-39","LT-40","LT-45")

female<-c("LT-17","LT-32","LT-35","LT-38","LT-41","LT-43",
          "LT-44","LT-46","LT-47")

unknown<-c("LT-10", "LT-11", "LT-12", "LT-13", "LT-14",
           "LT-15", "LT-16", "LT-18", "LT-19", "LT-2",
           "LT-20", "LT-21", "LT-24", "LT-25", "LT-26",
           "LT-27", "LT-28", "LT-29", "LT-30", "LT-4",
           "LT-5", "LT-7", "LT-9", "LT-23")

lt<-as.factor(c(male,female,unknown))

#assign tags to surgery groups
group2012<-paste0("LT-",c(1:30))
group2013<-paste0("LT-",c(31:44))
group2014<-paste0("LT-",c(45:47))

#assign tags with all data available for a given year
good2012<-c("LT-1","LT-2","LT-3","LT-4",
"LT-5","LT-7","LT-8","LT-9","LT-10","LT-11",
"LT-12","LT-13","LT-14","LT-15","LT-16","LT-17",
"LT-18","LT-19","LT-20","LT-21","LT-22","LT-24",
"LT-25","LT-26","LT-27","LT-28","LT-29","LT-30")

good2013<-c("LT-1","LT-2","LT-3","LT-5","LT-7",
"LT-8","LT-9","LT-13","LT-14","LT-15","LT-17",
"LT-18","LT-19","LT-21","LT-22","LT-26","LT-28",
"LT-29","LT-30","LT-31","LT-32","LT-33","LT-34",
"LT-35","LT-37","LT-40","LT-41","LT-43","LT-44")

good2014<-c("LT-1","LT-2","LT-3","LT-8","LT-9","LT-13",
"LT-14","LT-15","LT-17","LT-18","LT-19","LT-22",
"LT-26","LT-29","LT-30","LT-45","LT-46","LT-47")

  
##extract year,Month,day, hour

fish$YEAR<-year(fish$DATETIME)
fish$MONTH<-month(fish$DATETIME)
fish$DAY<-day(fish$DATETIME)
fish$HOUR<-hour(fish$DATETIME)

#assign sex M = male, F = female, U = unknown
fish$SEX<-"U"
fish$SEX[fish$TRANSMITTER %in% male]<-"M"
fish$SEX[fish$TRANSMITTER %in% female]<-"F"

#shore distance
#check if it exists
if(!exists("shore.dist")){
  #if not load shore dist asc raster
  require(SDMTools)
  shore.dist<-read.asc("~/Data/GIS/shoredist2.asc")
  #assign raster class object
  shore.dist<-raster.from.asc(shore.dist)
  #assign projection
  projection(shore.dist)<- "+proj=utm +zone=11 +ellps=WGS84"
  detach("package:SDMTools", unload=TRUE)
}


#extract distance from shore for each point
fish$SHOREDIST<-extract(shore.dist, fish[,c(1,2)])
fish$SHOREDIST[is.na(fish$SHOREDIST)]<-0

#extract bottom depth at fish location
fish$BOTTOMDEPTH<-(extract(r,fish[c("UTM.X","UTM.Y")]))

#calculate distance between bottom depth and fish depth
fish$OFFBOTTOM<-fish$BOTTOMDEPTH-fish$DEPTH

####NOT RUN
#####determine when 15 deg is less than 4m

#create temp and oxygen data frames
temp<-read.csv("~/Documents/Masters Thesis/Data/Temp/Alexie mean daily temp 2013.csv")
source('~/Documents/R-Directory/Masters Thesis/R-Functions/Functions/tempDepth.R')
require(reshape2)
temp<-melt(temp,id.var = "Depth..m.")

colnames(temp)<-c("Depth","Date","Temp")

#assign Date as date class (chron package)
temp$Date<-as.character(temp$Date)
temp$Date<-substr(temp$Date,2,11)
temp$Date<-gsub(".","/",temp$Date, fixed=TRUE)
temp$Date<-as.POSIXct(temp$Date, format = "%d / %m / %Y", tz="MST")

#remove NA's
temp<-temp[!is.na(temp$Date),]



temp2<-read.csv("~/Documents/Masters Thesis/Data/Temp/alexie_sept2013-may2014_temp.csv")

temp2$Date<-as.POSIXct(temp2$Date, format="%d / %m / %Y", tz="MST")
temp2<-temp2[,c("Depth","Date","Temp")]

temp<-rbind(temp,temp2)

save(temp,file="temp.RData")

depth15<-tempDepth(15,data=temp)

depth15$Depth[depth15$Depth==-1]<-0
#extend depth15 to october 15

#dates<-seq.POSIXt(from=as.POSIXct("2013-09-28",tz="MST"),
 #          to=as.POSIXct("2013-10-15",tz="MST"),
  #         by="day")
#tempext<-data.frame(Depth=rep(0,length(dates)),
 #                   Date=dates)
#depth15<-rbind(depth15,tempext)


depth15[depth15$Depth<=4,]

#spring defined as days before 15deg C thermocline falls deeper than 4m
#spawn is when  15deg C thermocline rises greater than 4m

######End of NOT RUN

#Extract fish positions during spawn ~ September 
fishspawn<-fish[fish$DATETIME>=as.POSIXct("2012-09-01",tz="MST") & 
                  fish$DATETIME<as.POSIXct("2012-10-01",tz="MST") |
                  fish$DATETIME>=as.POSIXct("2013-09-01",tz="MST") & 
                  fish$DATETIME<as.POSIXct("2013-10-01",tz="MST") |
                  fish$DATETIME>=as.POSIXct("2014-09-01",tz="MST") & 
                  fish$DATETIME<as.POSIXct("2014-10-01",tz="MST")
                ,]


clust.param<-data.frame(Year=NULL,
                        ID=NULL,
                        MinPts=NULL,
                        EpsSpace=NULL,
                        EpsTime=NULL)
year<-unique(fishspawn$YEAR)

#Loop for each year
for(i in 1:length(year)){
  #extract year specific data
  fishy<-fishspawn[fishspawn$YEAR==year[i],]
  #extract year specific tags
  if(year[i]==2012){
    tag<-good2012
  }else if (year[i]==2013){
    tag<-good2013
  }else
    tag<-good2014
  
  #Loop for each tag within a year
  for(ii in 1:length(tag)){
    #load fish
    fishi<-fishy[fishy$TRANSMITTER==tag[ii],]
    #order dataframe by increasing date
    fishi<-fishi[order(fishi$DATETIME),]
    #pull out required data x,y and time
    data<-data.frame(x=fishi$UTM.X,
                     y=fishi$UTM.Y,
                     t=fishi$DATETIME)
    
    #a good estimate for k or min points is the natural log of the number of data points
    k=round(log(dim(data)[1]))
    #create distance matrix for k nearest neighbours
    d<-get.knn(data[,c(1,2)],k)$nn.dist
    #take k-1 nearest neighbours for distance
    d<-d[,k-1]
    #order k-1 nearest neighbours by increasing order
    d<-d[order(d)]
    #create a time distance matrix for k nearest neighbours
    t<-get.knn(data[,3],k)$nn.dist
    #take k-1 nearest neighbours for time
    t<-t[,k-1]
    #order k-1 nearest neighbours by increasing order
    t<-t[order(t)]
    
    #plot k-dist for distance
    plot(d,type="l",main=paste(tag[ii],"K-Dist Space", sep=" "),xlab="Points sorted by distance",ylab=paste(k,"-NN distance",sep=""))
    #determine eps_space cuttoff
    hd<-round(mean(d)+2*sd(d))
    #draw cutoff line
    abline(h=hd,lty=2)
    
    #plot k-dist for time
    plot(t,type="l",main=paste(tag[ii],"K-Dist Time", sep=" "),xlab="Points sorted by distance",ylab=paste(k,"-NN distance",sep=""))
    #determine eps_time cuttoff
    ht<-round(mean(t)+2*sd(t))
    #draw cutoff line
    abline(h=ht,lty=2)
    
    clust.p<-data.frame(Year=year[i],
                        ID=tag[ii],
                        MinPts=k,
                        EpsSpace=hd,
                        EpsTime=ht)
    
    clust.param<-rbind(clust.param,clust.p)
  }
}
rm(fishy)
clust.param

                
#for each parameter a take the rounded mean of from the sampled animals for each year

#all
round(summary(clust.param$MinPts))
round(summary(clust.param$EpsSpace))
round(summary(clust.param$EpsTime))

#2012
min_pts_2012<-round(summary(clust.param$MinPts[clust.param$Year==2012]))
eps_space_2012<-round(summary(clust.param$EpsSpace[clust.param$Year==2012]))
eps_time_2012<-round(summary(clust.param$EpsTime[clust.param$Year==2012]))

#2013
min_pts_2013<-round(summary(clust.param$MinPts[clust.param$Year==2013]))
eps_space_2013<-round(summary(clust.param$EpsSpace[clust.param$Year==2013]))
eps_time_2013<-round(summary(clust.param$EpsTime[clust.param$Year==2013]))

#2013 - group 2012
min_pts_2013_g2012<-round(summary(clust.param$MinPts[clust.param$Year==2013 & clust.param$ID %in% group2012]))
eps_space_2013_g2012<-round(summary(clust.param$EpsSpace[clust.param$Year==2013 & clust.param$ID %in% group2012]))
eps_time_2013_g2012<-round(summary(clust.param$EpsTime[clust.param$Year==2013 & clust.param$ID %in% group2012]))

#2013 - group 2013
min_pts_2013_g2013<-round(summary(clust.param$MinPts[clust.param$Year==2013 & clust.param$ID %in% group2013]))
eps_space_2013_g2013<-round(summary(clust.param$EpsSpace[clust.param$Year==2013 & clust.param$ID %in% group2013]))
eps_time_2013_g2013<-round(summary(clust.param$EpsTime[clust.param$Year==2013 & clust.param$ID %in% group2013]))



#2014
min_pts_2014<-round(summary(clust.param$MinPts[clust.param$Year==2014]))
eps_space_2014<-round(summary(clust.param$EpsSpace[clust.param$Year==2014]))
eps_time_2014<-round(summary(clust.param$EpsTime[clust.param$Year==2014]))

#2014 - group 2012
min_pts_2014_g2012<-round(summary(clust.param$MinPts[clust.param$Year==2014 & clust.param$ID %in% group2012]))
eps_space_2014_g2012<-round(summary(clust.param$EpsSpace[clust.param$Year==2014 & clust.param$ID %in% group2012]))
eps_time_2014_g2012<-round(summary(clust.param$EpsTime[clust.param$Year==2014 & clust.param$ID %in% group2012]))

#2014 - group 2014
min_pts_2014_g2014<-round(summary(clust.param$MinPts[clust.param$Year==2014 & clust.param$ID %in% group2014]))
eps_space_2014_g2014<-round(summary(clust.param$EpsSpace[clust.param$Year==2014 & clust.param$ID %in% group2014]))
eps_time_2014_g2014<-round(summary(clust.param$EpsTime[clust.param$Year==2014 & clust.param$ID %in% group2014]))




#summarize all cluster parameters
cluster<-data.frame(Year=sort(rep(year,4)),
                    Group=rep(c("All","Group 2012","Group 2013","Group 2014"),3),
                    
                    Min_Pts_Median=c(min_pts_2012[["Median"]],min_pts_2012[["Median"]],"","",
                                     min_pts_2013[["Median"]],min_pts_2013_g2012[["Median"]],
                                     min_pts_2013_g2013[["Median"]],"",
                                     min_pts_2014[["Median"]],min_pts_2014_g2012[["Median"]],"",
                                     min_pts_2014_g2014[["Median"]]),
                    
                    Min_Pts_Mean=c(min_pts_2012[["Mean"]],min_pts_2012[["Mean"]],"","",
                                     min_pts_2013[["Mean"]],min_pts_2013_g2012[["Mean"]],
                                     min_pts_2013_g2013[["Mean"]],"",
                                     min_pts_2014[["Mean"]],min_pts_2014_g2012[["Mean"]],"",
                                     min_pts_2014_g2014[["Mean"]]),
                    
                    Min_Pts_Min=c(min_pts_2012[["Min."]],min_pts_2012[["Min."]],"","",
                                  min_pts_2013[["Min."]],min_pts_2013_g2012[["Min."]],min_pts_2013_g2013[["Min."]],"",
                                  min_pts_2014[["Min."]],min_pts_2014_g2012[["Min."]],"",min_pts_2014_g2014[["Min."]]),
                    
                    Min_Pts_Max=c(min_pts_2012[["Max."]],min_pts_2012[["Max."]],"","",
                                  min_pts_2013[["Max."]],min_pts_2013_g2012[["Max."]],min_pts_2013_g2013[["Max."]],"",
                                  min_pts_2014[["Max."]],min_pts_2014_g2012[["Max."]],"",min_pts_2014_g2014[["Max."]]),
                    
                    Eps_Space_Median=c(eps_space_2012[["Median"]],eps_space_2012[["Median"]],"","",
                                       eps_space_2013[["Median"]],eps_space_2013_g2012[["Median"]],
                                       eps_space_2013_g2013[["Median"]],"",
                                       eps_space_2014[["Median"]],eps_space_2014_g2012[["Median"]],"",
                                      eps_space_2014_g2014[["Median"]]),
                    
                    Eps_Space_Mean=c(eps_space_2012[["Mean"]],eps_space_2012[["Mean"]],"","",
                                       eps_space_2013[["Mean"]],eps_space_2013_g2012[["Mean"]],
                                       eps_space_2013_g2013[["Mean"]],"",
                                       eps_space_2014[["Mean"]],eps_space_2014_g2012[["Mean"]],"",
                                       eps_space_2014_g2014[["Mean"]]),
                    
                    Eps_Space_Min=c(eps_space_2012[["Min."]],eps_space_2012[["Min."]],"","",
                                    eps_space_2013[["Min."]],eps_space_2013_g2012[["Min."]],
                                    eps_space_2013_g2013[["Min."]],"",
                                    eps_space_2014[["Min."]],eps_space_2014_g2012[["Min."]],"",
                                    eps_space_2014_g2014[["Min."]]),
                    
                    Eps_Space_Max=c(eps_space_2012[["Max."]],eps_space_2012[["Max."]],"","",
                                    eps_space_2013[["Max."]],eps_space_2013_g2012[["Max."]],
                                    eps_space_2013_g2013[["Max."]],"",
                                    eps_space_2014[["Max."]],eps_space_2014_g2012[["Max."]],"",
                                    eps_space_2014_g2014[["Max."]]),
                    
                    Eps_Time_Median=c(eps_time_2012[["Median"]],eps_time_2012[["Median"]],"","",
                                      eps_time_2013[["Median"]],eps_time_2013_g2012[["Median"]],
                                      eps_time_2013_g2013[["Median"]],"",
                                      eps_time_2014[["Median"]],eps_time_2014_g2012[["Median"]],"",
                                      eps_time_2014_g2014[["Median"]]),
                    
                    Eps_Time_Mean=c(eps_time_2012[["Mean"]],eps_time_2012[["Mean"]],"","",
                                      eps_time_2013[["Mean"]],eps_time_2013_g2012[["Mean"]],
                                      eps_time_2013_g2013[["Mean"]],"",
                                      eps_time_2014[["Mean"]],eps_time_2014_g2012[["Mean"]],"",
                                      eps_time_2014_g2014[["Mean"]]),
                    
                    Eps_Time_Min=c(eps_time_2012[["Min."]],eps_time_2012[["Min."]],"","",
                                   eps_time_2013[["Min."]],eps_time_2013_g2012[["Min."]],
                                   eps_time_2013_g2013[["Min."]],"",
                                   eps_time_2014[["Min."]],eps_time_2014_g2012[["Min."]],"",
                                   eps_time_2014_g2014[["Min."]]),
                    
                    Eps_Time_Max=c(eps_time_2012[["Max."]],eps_time_2012[["Max."]],"","",
                                   eps_time_2013[["Max."]],eps_time_2013_g2012[["Max."]],
                                   eps_time_2013_g2013[["Max."]],"",
                                   eps_time_2014[["Max."]],eps_time_2014_g2012[["Max."]],"",
                                   eps_time_2014_g2014["Max."])
)

#define cluster parameters taken as mean of new groups 2012 group and 2013 group
#2012 group data deteriorates over time and 2014 tags have different settings than 2012 and 2013 groups
min_pts<-mean(min_pts_2012[["Mean"]],min_pts_2013_g2013[["Mean"]])
eps_space<-mean(eps_space_2012[["Mean"]],eps_space_2013_g2013[["Mean"]])
eps_time<-mean(eps_time_2012[["Mean"]],eps_time_2013_g2013[["Mean"]])



 
#CREATE CLUSTERS
year<-unique(fish$YEAR)

#Loop for each year
for(i in 1:length(year)){
  
  #extract year specific data
  fishy<-fish[fish$YEAR==year[i],]
 
   #extract year specific tags
  if(year[i]==2012){
    tag<-good2012
  }else if (year[i]==2013){
    tag<-good2013
  }else
    tag<-good2014
  #loop for each good tag in each year
  for(ii in 1:length(tag)){
    #extract data for specific tag within given year
    fishi<-fishy[fishy$TRANSMITTER==tag[ii],]
    #order data by timestamp
    fishi<-fishi[order(fishi$DATETIME),]
    #extract variables for stdbscan function
    x=fishi$UTM.X
    y=fishi$UTM.Y
    t=fishi$DATETIME
    #run stdbscan function
    clust<-stdbscan(lon=x, lat=y, time=t, eps_space=eps_space, eps_time=eps_time, deltaE=2,SD=TRUE,MinPts=min_pts)
    
    #print status of running process: year  - progress % within year - tag ID 
    print(year[i])
    #print(clust)
    print(paste0("Progress = ", round((ii/length(tag))*100),"%"))
    print(tag[ii])
    #assign cluster id 
    fishi$CLUSTER<-clust$cluster
    assign(paste("LT",substr(tag[ii],start=4,stop=5),"y",year[i],sep=""),fishi)
    assign(paste("clusterLT",substr(tag[ii],start=4,stop=5),"y",year[i],sep=""),clust)
  }
}



fishlist = lapply(ls()[grep("(?<!cluster)LT",ls(),perl=TRUE)], get)
clustlist = lapply(ls()[grep("clusterLT",ls())], get)

#remove individual data frames that are not needed
rm(fishi)
rm(list=ls()[grep("(?<!cluster)LT",ls(),perl=TRUE)])
rm(list=ls()[grep("clusterLT",ls())])

#need spawn_pts and egg_pts from All Spawn Sites.R
#egg_pts<-read.csv("~/Data/GIS/eggsites.csv")

#Load spawning site locations (both egg present and unknown)
spawn_pts<-read.csv("~/Data/GIS/spawnsites.csv")
#create spatial points and project them
coordinates(spawn_pts)<-~UTM.X+UTM.Y
proj4string(spawn_pts)<-CRS("+proj=utm +zone=11 +ellps=WGS84")

#Create 20m buffer around points
spawnpoly<-gBuffer(spawn_pts,byid=TRUE,width=20)



for(i in 1:length(fishlist)){
  spawnclust=NULL
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))
  plot(shore_outline)
  year<-fishi$YEAR[[1]]
  tag<-fishi$TRANSMITTER[1]
  title(paste(year,tag,sep=" "))
  
  lines(fishi$UTM.X,fishi$UTM.Y,col=rgb(0,0,0,0.2),pch=19,cex=0.2)
  points(xy,col=rgb(0,0,0,0.2),pch=19,cex=0.02)
  if (max(fishi$CLUSTER)>0){
    for(ii in 1:max(fishi$CLUSTER)){
      new.xy<-xy[fishi$CLUSTER==ii]
      colour<-"red"
      cp=mcp(new.xy,percent=100)
      
      plot(cp,add=TRUE,col=colour)
      
      spawnclust2<-data.frame(TRANSMITTER=tag,
                              YEAR=year,
                              CLUSTER=ii,
                              X=mean(cp@bbox[1,c(1,2)]),
                              Y=mean(cp@bbox[2,c(1,2)]),
                              AREA=gArea(cp),
                              SPAWNSITE=ifelse(any(gOverlaps(cp,spawnpoly, byid=TRUE)),1,0),
                              EGGS=ifelse(any(gOverlaps(cp,spawnpoly[spawnpoly$EGGS==1,], byid=TRUE)),1,0))
      
      spawnclust<-rbind(spawnclust,spawnclust2)
      
    }
    assign(paste("spawnclustLT",substr(tag,start=4,stop=5),"y",year,sep=""),spawnclust)
  }

}

#Create a cluster list NOTE 72 out of 75 create a cluster 
clustlist<-lapply(ls()[grep("spawnclustLT",ls())], get)
#remove individual data.frames
rm(list=ls()[grep("spawnclust",ls())])


#Spawn cluster filter
#create new dataframe for lake trout cluster data
ltspawnclust<-NULL

#assign "year tag" name to each list object 
fishlist<-setNames(fishlist,lapply(fishlist, function(x) paste0(x$YEAR[1],x$TRANSMITTER[1])))

#vector of accel "year tag" names
acceltag<-paste0(rep(2013,14),"LT-",c(31:44))
#find accel tags in fishlist
acceltag<-names(fishlist)[names(fishlist) %in% acceltag]

#determine accel cutoff
fishbind<-do.call(rbind,fishlist[acceltag])
fishbind_all<-do.call(rbind,fishlist)
fishbind_sept<-fishbind_all[fishbind_all$MONTH==9,]

#mean ACCEL+2 sd =accel of 1.3
accel_cut<-mean(fishbind$ACCEL[fishbind$DATETIME>="2013-08-30" 
                               & fishbind$DATETIME<"2013-10-01" 
                               & fishbind$CLUSTER!=0], na.rm=TRUE) +
                2 * (sd(fishbind$ACCEL[fishbind$DATETIME>="2013-08-30" 
                                       & fishbind$DATETIME<"2013-10-01" 
                                       & fishbind$CLUSTER!=0],na.rm=TRUE))


for(i in 1:length(clustlist)){
  
  #extract clustlist object
  spawnclust=clustlist[[i]]
  
  #determine year and tag ID
  tag=spawnclust$TRANSMITTER[1]
  year=spawnclust$YEAR[1]
  yeartag=paste0(year,tag)
  
  #Extract matching position data
  fishi=fishlist[[yeartag]]
  
  if(length(spawnclust$CLUSTER)>0){
    #extract fish posiiton data inside a cluster
    fishi<-fishi[fishi$CLUSTER %in% spawnclust$CLUSTER,]
    
    #create cluster summary
    for(ii in 1: length(spawnclust$CLUSTER)){
      fishi2<-fishi[fishi$CLUSTER==spawnclust$CLUSTER[ii],]
      ltspawnclust2<-data.frame(TRANSMITTER=fishi2$TRANSMITTER[1],
                                YEAR=fishi2$YEAR[1],
                                CLUSTER=spawnclust$CLUSTER[ii],
                                #centre position of cluster
                                CLUST.X=spawnclust$X[ii],
                                CLUST.Y=spawnclust$Y[ii],
                                #area of cluster
                                AREA=spawnclust$AREA[ii],
                                #binary if it overlaps a possible spawning site; 1 = true
                                SPAWNSITE=spawnclust$SPAWNSITE[ii],
                                #binary if it overlaps a spawn site with eggs present; 1 = true
                                EGGS=spawnclust$EGGS[ii],
                                #start of cluster
                                STARTTIME=min(fishi2$DATETIME),
                                #duration of cluster
                                TOTALTIME=as.numeric(difftime(max(fishi2$DATETIME),min(fishi2$DATETIME),units="mins")),
                                #mean depth of fish in cluster
                                DEPTH=mean(fishi2$DEPTH,na.rm=TRUE),
                                #min depth of fish in cluster
                                MINDEPTH=min(fishi2$DEPTH,na.rm=TRUE),
                                #max depth of fish in cluster
                                MAXDEPTH=max(fishi2$DEPTH,na.rm=TRUE),
                                #mean lake depth at fish posiitons
                                BOTTOMDEPTH=mean(fishi2$BOTTOMDEPTH,na.rm=TRUE),
                                #percent of positions at lake depths less than 4 meters
                                PERCUNDER4M=length(fishi2$BOTTOMDEPTH[fishi2$BOTTOMDEPTH<4])/length(fishi2$BOTTOMDEPTH),
                                #max lake depth at fish positions
                                MAXBOTTOMDEPTH=max(fishi2$BOTTOMDEPTH,na.rm=TRUE),
                                #min lake depth at fish positions
                                MINBOTTOMDEPTH=min(fishi2$BOTTOMDEPTH,na.rm=TRUE),
                                #mean distance of positions from shore
                                SHOREDIST=mean(fishi2$SHOREDIST,na.rm=TRUE),
                                #percent positions less than 50 meters from shore
                                PERCUNDER50M=length(fishi2$SHOREDIST[fishi2$SHOREDIST<50])/length(fishi2$SHOREDIST),
                                #mean acceleration of fish in cluster
                                ACCEL=mean(fishi2$ACCEL,na.rm=TRUE),
                                #max acceleration of fish in cluster
                                MAXACCEL=max(fishi2$ACCEL,na.rm=TRUE),
                                #percent over acceleration cutoff of 1.3 for positions in cluster
                                PERCOVERACCELCUT=length(fishi2$ACCEL[fishi2$ACCEL>accel_cut & !is.na(fishi2$ACCEL)])/length(fishi2$ACCEL[!is.na(fishi2$ACCEL)])
      )
      ltspawnclust<-rbind(ltspawnclust,ltspawnclust2)
    }
  }
}

#assign sex
#U for unknown
ltspawnclust$SEX<-"U"
#M for male
ltspawnclust$SEX[ltspawnclust$TRANSMITTER %in% male]<-"M"
#F for female
ltspawnclust$SEX[ltspawnclust$TRANSMITTER %in% female]<-"F"

plot(shore_outline,main="All Clusters")
points(CLUST.Y~CLUST.X,data=ltspawnclust,pch=19,cex=0.3,col=rgb(0,0,1,0.1))
#area_cut=mean(ltspawnclust$AREA)+2*sd(ltspawnclust$AREA)
#data1=ltspawnclust[ltspawnclust$AREA<area_cut,]


#plot(shore_outline,main="Potenital Spawn Clusters")
#select data with >50 percent positions under 50 m distance to shore 
data1=ltspawnclust[ltspawnclust$PERCUNDER50M>.50,]
#points(CLUST.Y~CLUST.X,data=data1,pch=19,col=rgb(1,0,0,0.3))

#select data with >50 percent positions at a lake bottom depth of less than 4 m
data2=data1[data1$PERCUNDER4M>.50,]
#points(CLUST.Y~CLUST.X,data=data2,pch=19,cex=0.5,col=rgb(0,0,1,1))

#create a new dataframe
ltspawnclust2<-data2

#Determine unique spawning sites ie. greater than 70 meters apart
ltspawnclust2$SITE<-NA

#all unique sites
for(i in 1:dim(ltspawnclust2)[1]){
  if(is.na(ltspawnclust2$SITE[i])){
    d<-pointDistance(ltspawnclust2[i,c("CLUST.X","CLUST.Y")],ltspawnclust2[,c("CLUST.X","CLUST.Y")],allpairs=TRUE,longlat=FALSE)
    #the centre of clusters must be 70 m apart to be concidered "unique sites"
    ltspawnclust2$SITE[d<70]<-i
  }
}


ltspawnclust<-ltspawnclust2
rm(ltspawnclust2)

ltspawnclust$SPAWN<-0
ltspawnclust$SPAWN[ltspawnclust$MAXACCEL > accel_cut]<-1

#add date time columns
ltspawnclust$HOUR<-hour(ltspawnclust$STARTTIME)+minute(ltspawnclust$STARTTIME)/60+second(ltspawnclust$STARTTIME)/3600
ltspawnclust$DAY<-day(ltspawnclust$STARTTIME)
ltspawnclust$MONTH<-month(ltspawnclust$STARTTIME)
#reorder rownames
rownames(ltspawnclust)<-seq(1,dim(ltspawnclust)[1],by=1)

#plots 6.15 x 4.76
plot(shore_outline)
points(CLUST.Y~CLUST.X,data=ltspawnclust)
points(CLUST.Y~CLUST.X,data=ltspawnclust[ltspawnclust$SPAWN==1,],pch=19,col=rgb(1,0,0,1))
points(CLUST.Y~CLUST.X,data=ltspawnclust[ltspawnclust$SPAWNSITE==1,],pch=19,col=rgb(0,1,0,1))
points(CLUST.Y~CLUST.X,data=ltspawnclust[ltspawnclust$EGGS==1,],pch=19,col=rgb(0,0,1,1))


clust<-unlist(dimnames(table(ltspawnclust$SITE)[table(ltspawnclust$SITE)]))

as.numeric(table(ltspawnclust$SITE))
plot(table(ltspawnclust$SITE),type="h")

plot(shore_outline)
points(CLUST.Y~CLUST.X,data=ltspawnclust[ltspawnclust$SITE %in% clust,],col=YEAR)


plot(shore_outline)
points(CLUST.Y~CLUST.X,data=ltspawnclust[ltspawnclust$YEAR==2014,],pch=19,col=rgb(1,0,0,0.1))



ltspawnclust$DATE<-format(ltspawnclust$STARTTIME,"%Y-%m-%d")

#all clusters
spawn.num=NULL
spawn.num$DATE<-as.POSIXct(names(table(ltspawnclust$DATE)),format="%Y-%m-%d", tz="MST")
spawn.num$YEAR<-year(spawn.num$DATE)
spawn.num$MONTH<-month(spawn.num$DATE)
spawn.num$FREQ<-as.numeric(table(ltspawnclust$DATE))
spawn.num$FREQ12<-as.numeric(table(ltspawnclust$DATE[ltspawnclust$TRANSMITTER %in% group2012]))
spawn.num$FREQ13<-as.numeric(table(ltspawnclust$DATE[ltspawnclust$TRANSMITTER %in% group2013]))
spawn.num$FREQ14<-as.numeric(table(ltspawnclust$DATE[ltspawnclust$TRANSMITTER %in% group2014]))
spawn.num<-as.data.frame(spawn.num)

length(unique(ltspawnclust$TRANSMITTER[ltspawnclust$TRANSMITTER %in% group2012]))
length(unique(ltspawnclust$TRANSMITTER[ltspawnclust$TRANSMITTER %in% group2013]))
length(unique(ltspawnclust$TRANSMITTER[ltspawnclust$TRANSMITTER %in% group2014]))


#UTM 
data<-ltspawnclust[ltspawnclust$MONTH==9,c(2,4,5,24)]

library(dplyr)
df2 <- data %>% 
  group_by(SITE, YEAR) %>%
  summarise(FREQ=length(SITE))
df3 <- data %>%
  group_by(SITE) %>%
  summarise(CLUST.X=mean(CLUST.X),
            CLUST.Y=mean(CLUST.Y))

first.clust <- df2 %>%
  left_join(df3)




library(mapplots)

data<-first.clust
xyz<-make.xyz(x=data$CLUST.X,y=data$CLUST.Y,z=data$FREQ,group=data$YEAR)

plot(shore_outline)
  
draw.pie(z=xyz$z, x=xyz$x, y=xyz$y, 
        radius=100, col=c(alpha("blue", 0.8), alpha("red", 0.8), alpha("yellow", 0.8)),
        border=c(alpha("blue", 0.8), alpha("red", 0.8), alpha("yellow", 0.8))) 

text(650800,6954500, "Year",font=2) 
legend.pie(650800,6954300, labels=c("2012","2013","2014"), radius = 100,
           bty="n", col=c(alpha("blue", 0.8), alpha("red", 0.8), alpha("yellow", 0.8)),
           border=c(alpha("blue", 0.8), alpha("red", 0.8), alpha("yellow", 0.8)),
           label.dist=1.7)

text(650800,6954000, "Frequency",font=2) 
legend.bubble(650800,6953800, z=92, round=1, n=1, maxradius=100, bty="n",txt.col="transparent")
text(650950,6953800, "92", cex=0.8) 
legend.bubble(650800,6953625, z=46, round=1, n=1, maxradius=50, bty="n",txt.col="transparent")
text(650900,6953625, "46", cex=0.7) 
legend.bubble(650800,6953525, z=23, round=1, n=1, maxradius=25, bty="n",txt.col="transparent")
text(650875,6953525,"23", cex=0.6) 
legend.bubble(650800,6953462, z=12, round=1, n=1, maxradius=13, bty="n",txt.col="transparent")
text(650862,6953462, "12", cex=0.5)               


#determine max frequency
max.all<-max(first.clust$FREQ)

#2012
data<-filter(first.clust,YEAR==2012)
#determine scaling factor
s=max(data$FREQ)/max.all

plot(shore_outline, main=data$YEAR[1], col="grey90",border="grey70")

draw.bubble(z=data$FREQ, x=data$CLUST.X, y=data$CLUST.Y, 
         maxradius=200*s, col=alpha("lightblue", 0.3),pch=19)
draw.bubble(z=data$FREQ, x=data$CLUST.X, y=data$CLUST.Y, 
            maxradius=200*s, col=alpha("blue",0.3))


legend.bubble(650800,6954000, z=max.all, round=1, n=3, maxradius=200, bty="n",pch=19,
              col=alpha("lightblue",0.3), txt.cex = 0.7, txt.col="transparent")
legend.bubble(650800,6954000, z=max.all, round=0, n=3, maxradius=200, bty="n",
              col=alpha("blue",0.3),txt.cex = 0.7, txt.col = alpha("blue",0.3))
#2013
data<-filter(first.clust,YEAR==2013)
#determine scaling factor
s=max(data$FREQ)/max.all
plot(shore_outline, main=data$YEAR[1], col="grey90",border="grey70")

draw.bubble(z=data$FREQ, x=data$CLUST.X, y=data$CLUST.Y, 
            maxradius=200*s, col=alpha("lightblue", 0.3),pch=19)
draw.bubble(z=data$FREQ, x=data$CLUST.X, y=data$CLUST.Y, 
            maxradius=200*s, col=alpha("blue",0.3))


legend.bubble(650800,6954000, z=max.all, round=1, n=3, maxradius=200, bty="n",pch=19,
              col=alpha("lightblue",0.3), txt.cex = 0.7, txt.col="transparent")
legend.bubble(650800,6954000, z=max.all, round=0, n=3, maxradius=200, bty="n",
              col=alpha("blue",0.3),txt.cex = 0.7, txt.col = alpha("blue",0.3))
#2014
data<-filter(first.clust,YEAR==2014)
#determine scaling factor
s=max(data$FREQ)/max.all

plot(shore_outline, main=data$YEAR[1], col="grey90",border="grey70")

draw.bubble(z=data$FREQ, x=data$CLUST.X, y=data$CLUST.Y, 
            maxradius=200*s, col=alpha("lightblue", 0.3),pch=19)
draw.bubble(z=data$FREQ, x=data$CLUST.X, y=data$CLUST.Y, 
            maxradius=200*s, col=alpha("blue",0.3))


legend.bubble(650800,6954000, z=max.all, round=1, n=3, maxradius=200, bty="n",pch=19,
              col=alpha("lightblue",0.3), txt.cex = 0.7, txt.col="transparent")
legend.bubble(650800,6954000, z=max.all, round=0, n=3, maxradius=200, bty="n",
              col=alpha("blue",0.3),txt.cex = 0.7, txt.col = alpha("blue",0.3))

data<-table(ltspawnclust$YEAR,ltspawnclust$SITE)
data<-data[,order(colSums(data))]



barplot(data,horiz=TRUE,col=c(alpha("blue",.10),alpha("blue",.50),alpha("blue",1)),
        border=alpha("blue",0.3))



peakdates<-function(date,freq){
  
  p<-data.frame(DATE=NULL,
                COUNT=NULL,
                DAYS=NULL)
  
  data<-data.frame(DATE=date,
                    FREQ=freq)
  
  for(i in 1:dim(data)[1]){
    data2<-data[i:dim(data)[1],]
    c50<-sum(data$FREQ)*0.5
    count<-sum(data2$FREQ[which(!cumsum(data2$FREQ)>c50)])
    days<-length(which(!cumsum(data2$FREQ)>c50))
    end<-dim(data2)[1]
    
    if(end>days){
      if((count+data2$FREQ[days+1])>c50){
        days=days
      }else{
        days=days=0
      }
      
    }else{
      days=0
    }
    p<-rbind(p,data.frame(DATE=data$DATE[i],COUNT=count,DAYS=days))
  }
  peak<-data.frame(START=min(p$DATE[p$DAYS==min(p$DAYS[p$DAYS>0])]),
                   END=min(p$DATE[p$DAYS==min(p$DAYS[p$DAYS>0])])+min(p$DAYS[p$DAYS>0])*24*60*60)
  print(peak)
}

peakdates.12<-peakdates(date=spawn.num.12$DATE,freq=spawn.num.12$FREQ)
peakdates.13<-peakdates(date=spawn.num.13$DATE,freq=spawn.num.13$FREQ)
peakdates.g13<-peakdates(date=spawn.num.13$DATE,freq=spawn.num.13$FREQ13)
peakdates.14<-peakdates(date=spawn.num.14$DATE,freq=spawn.num.14$FREQ)
peakdates.g14<-peakdates(date=spawn.num.14$DATE,freq=spawn.num.14$FREQ14)
op<-par()
par(mfrow=c(3,1))

#2012
plot(FREQ/24~DATE,data=spawn.num.12,ylim=c(0,3),
     xlim=as.POSIXct(c("2012-08-15","2012-10-15"),format="%Y-%m-%d", tz="MST"),type="h",col="grey",
     lwd=6,lend=2)
lines(c(peakdates.12$START-12*60*60,peakdates.12$START-12*60*60),c(-1,3.5))
lines(c(peakdates.12$END+12*60*60,peakdates.12$END+12*60*60),c(-1,3.5))

#2013 all
plot(FREQ~DATE,data=spawn.num.13,ylim=c(0,50),
     xlim=as.POSIXct(c("2013-08-15","2013-10-15"),format="%Y-%m-%d", tz="MST"), type="h",col="grey",
     lwd=6,lend=2)
lines(c(peakdates.13$START-12*60*60,peakdates.13$START-12*60*60),c(-1,50))
lines(c(peakdates.13$END+12*60*60,peakdates.13$END+12*60*60),c(-1,50))


#2014 all
plot(FREQ~DATE,data=spawn.num.14,ylim=c(0,10),
     xlim=as.POSIXct(c("2014-08-15","2014-10-15"),format="%Y-%m-%d", tz="MST"),type="h",col="grey",
     lwd=6,lend=2)
lines(c(peakdates.14$START-12*60*60,peakdates.14$START-12*60*60),c(-1,10))
lines(c(peakdates.14$END+12*60*60,peakdates.14$END+12*60*60),c(-1,10))

par(op)

#2013 group 2013
plot(FREQ13/9~DATE,data=spawn.num.13,ylim=c(0,3),
     xlim=as.POSIXct(c("2013-08-15","2013-10-15"),format="%Y-%m-%d", tz="MST"), type="h",col="grey",
     lwd=6,lend=2)
lines(c(peakdates.g13$START-12*60*60,peakdates.g13$START-12*60*60),c(-1,3.5))
lines(c(peakdates.g13$END+12*60*60,peakdates.g13$END+12*60*60),c(-1,3.5))

#2014 group 2014
plot(FREQ14/3~DATE,data=spawn.num.14,ylim=c(0,3),
     xlim=as.POSIXct(c("2014-08-15","2014-10-15"),format="%Y-%m-%d", tz="MST"),type="h",col="grey",
     lwd=6,lend=2)
lines(c(peakdates.g14$START-12*60*60,peakdates.g14$START-12*60*60),c(-1,3))
lines(c(peakdates.g14$END+12*60*60,peakdates.g14$END+12*60*60),c(-1,3))




#daily movements in September
#run through each list segement (tag year)

#add steplength, turn angle to fishlist objects
for(i in 1:length(fishlist)){
  data<-fishlist[[i]]
  x<-data$UTM.X
  y<-data$UTM.Y
  #add x and y data as real and imaginary numbers in the z vector
  z<-x+1i*y
  
  #determine step lengths
  dz<-diff(z)
  #steplength
  s<-Mod(dz)
  #start from 2nd position because you require 2 positions to get a step length
  data$STEPLENGTH[2:length(z)]<-s
  
  #absolute angles
  phi<-Arg(dz)
  #turning angles
  theta<-diff(phi)
  #start from position 3 because you require 3 positions to get a turning angle
  data$TURNANGLE[3:length(z)]<-theta
  data$DATE<-format(data$DATETIME,"%Y-%m-%d")
  fishlist[[i]]<-data
}

fishbind.all<-bind_rows(fishlist)
fishbind.all$YEAR<-as.factor(fishbind.all$YEAR)
fishbind.all$DATE<-as.factor(fishbind.all$DATE)

daily <- fishbind.all %>%
  group_by(TRANSMITTER,YEAR,DATE) %>%
  summarise(MONTH=MONTH[1],
            SEX=SEX[1],
            STEPLENGTH=sum(STEPLENGTH, na.rm=TRUE),
            PERSINDEX=mean(cos(TURNANGLE), na.rm=TRUE),
            ACCELMEAN=mean(ACCEL,na.rm=TRUE),
            ACCELMEDIAN=median(ACCEL,na.rm=TRUE),
            ACCELMAX=max(ACCEL,na.rm=TRUE),
            ACCELSUM=sum(ACCEL,na.rm=TRUE))
  
#remove download dates (partial data)

daily<-daily[daily$DATE!="2012-09-11" & daily$DATE!="2012-09-13",]
daily<-daily[daily$DATE!="2013-09-24" & daily$DATE!="2013-09-27",]

daily$ACCELMEAN[is.nan(daily$ACCELMEAN)]<-NA

spawntable<- daily %>%
  #select(MONTH==9) %>%
    group_by(YEAR,MONTH,SEX) %>%
      summarise(n = length(unique(TRANSMITTER)),
                STEPLENGTH.MEAN=mean(STEPLENGTH, na.rm=TRUE),
                STEPLENGTH.SD=sd(STEPLENGTH, na.rm=TRUE),
                PERSINDEX.MEAN=mean(PERSINDEX, na.rm=TRUE),
                PERSINDEX.SD=sd(PERSINDEX, na.rm=TRUE),
                ACCEL.MEAN=mean(ACCELMEAN,na.rm=TRUE),
                ACCEL.SD=sd(ACCELMEAN,na.rm=TRUE))
library(knitr)

kable(spawntable[spawntable$MONTH==9,])




#distance between spawning shoals
spawn<-NULL
data<-ltspawnclust[ltspawnclust$MONTH==9,]
years<-unique(data$YEAR)
#loop for each year
for(i in 1:length(years)){
  #extract year data
  data.year<-data[data$YEAR==years[i],]
  #extract unique tags
  tags<-unique(data.year$TRANSMITTER)
  #loop for each unique tag
  for(ii in 1:length(tags)){
    #extract tag specific data
    fishi<-data.year[data.year$TRANSMITTER==tags[ii],]
    fishi$DISTCLUST<-NA
    fishi$TIMECLUST<-NA
    
    if(dim(fishi)[1]>1){
      t<-fishi$STARTTIME
      x<-fishi$CLUST.X
      y<-fishi$CLUST.Y
      #add x and y data as real and imaginary numbers in the z vector
      z<-x+1i*y
      #determine step lengths
      fishi$DISTCLUST[2:length(z)]<-Mod(diff(z))
      fishi$TIMECLUST[2:length(z)]<-difftime(t[c(2:length(t))],t[c(1:length(t)-1)],units="hours")
    }
    
    fishi$DISTSPAWNCLUST=NA
    fishi$TIMESPAWNCLUST=NA 
    fishi$CLUSTPRIOR=NA
    fishi2<-fishi[fishi$SPAWN==1,]
    if(dim(fishi2)[1]>0){
      ## clusters prior to first high activity cluster
      fishi$CLUSTPRIOR<-match(1,fishi$SPAWN)-1
      if(dim(fishi2)[1]>1){
        t<-fishi2$STARTTIME
        x<-fishi2$CLUST.X
        y<-fishi2$CLUST.Y
        #add x and y data as real and imaginary numbers in the z vector
        z<-x+1i*y
        fishi[fishi$SPAWN==1,]$DISTSPAWNCLUST[2:length(z)]<-Mod(diff(z))
        fishi[fishi$SPAWN==1,]$TIMESPAWNCLUST[2:length(z)]<-difftime(t[c(2:length(t))],t[c(1:length(t)-1)],units="hours")
      }
    }
    fishi$DISTNOCLUST=NA
    fishi$TIMENOCLUST=NA 
    fishi2<-fishi[fishi$SPAWN==0,]
    if(dim(fishi2)[1]>1){
      t<-fishi2$STARTTIME
      x<-fishi2$CLUST.X
      y<-fishi2$CLUST.Y
      #add x and y data as real and imaginary numbers in the z vector
      z<-x+1i*y
      fishi[fishi$SPAWN==0,]$DISTNOCLUST[2:length(z)]<-Mod(diff(z))
      fishi[fishi$SPAWN==0,]$TIMENOCLUST[2:length(z)]<-difftime(t[c(2:length(t))],t[c(1:length(t)-1)],units="hours")
    }
    
    spawn<-rbind(spawn,fishi)
  }
}


#UTM 



require(adehabitatHR)
mcpclust=NULL
years<-unique(spawn$YEAR)

for(i in 1:length(years)){
  spawn.year<-spawn[spawn$YEAR==years[i],]
  #extract unique tags
  tags<-unique(spawn.year$TRANSMITTER)
  #loop for each unique tag
  for(ii in 1:length(tags)){
    #extract tag specific data
    fishi<-spawn.year[spawn.year$TRANSMITTER==tags[ii],]
    
    xy<-SpatialPoints(cbind(fishi$CLUST.X,fishi$CLUST.Y))
    cp100=NULL
    cp50=NULL
    cp100$area=NA
    cp50$area=NA
    if(length(xy)>4){
      cp100=adehabitatHR::mcp(xy,percent=100,unout="km2")
      cp50=adehabitatHR::mcp(xy,percent=50,unout="km2")
    }
    #spawn clusters
    fishi2<-fishi[fishi$SPAWN==1,]
    cps=NULL
    cps$area=NA
    if(dim(fishi2)[1]>4){
      new.xy<-SpatialPoints(cbind(fishi2$CLUST.X,fishi2$CLUST.Y))
      cps=adehabitatHR::mcp(new.xy,percent=100,unout="km2")
    }
    fishi2<-fishi[fishi$SPAWN==0,]
    cpn=NULL
    cpn$area=NA
    if(dim(fishi2)[1]>4){
      new.xy<-SpatialPoints(cbind(fishi2$CLUST.X,fishi2$CLUST.Y))
      cpn=adehabitatHR::mcp(new.xy,percent=100,unout="km2")
    }
    mcpclust1<-data.frame(TRANSMITTER=tags[ii],
                          SEX=fishi$SEX[1],
                          MCP100=cp100$area,
                          MCP50=cp50$area,
                          SPAWNMCP=cps$area,
                          NOMCP=cpn$area,
                          PERC=cps$area/cp100$area*100)
    mcpclust<-rbind(mcpclust,mcpclust1)
    assign(paste("mcp100lt",substr(tags[ii],start=4,stop=5),"y",years[i],sep=""),cp100)
    assign(paste("mcp50lt",substr(tags[ii],start=4,stop=5),"y",years[i],sep=""),cp50)
  }
}
 
require(adehabitatHR)
xy<-SpatialPoints(cbind(spawn$CLUST.X,spawn$CLUST.Y))
all_mcp<-adehabitatHR::mcp(xy,percent=100,unout="km2")

####overlap
?gOverlaps

save(ltspawnclust,first.clust,spawn.num,cluster,shore_outline,
     daily, spawntable, spawn, mcpclust, all_mcp, file="SpawnClust.RData")

