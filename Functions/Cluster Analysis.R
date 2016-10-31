
#plots 6.15 x 4.76
source('~/Google Drive/st-dbscan/st-dbscan function.R')
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

male<-c("LT-1","LT-3","LT-8","LT-22","LT-31","LT-33",
        "LT-34","LT-36","LT-37","LT-39","LT-40","LT-42",
        "LT-45")
female<-c("LT-17","LT-32","LT-35","LT-38","LT-41","LT-43",
          "LT-44","LT-46","LT-47")
unknown<-c("LT-10", "LT-11", "LT-12", "LT-13", "LT-14",
           "LT-15", "LT-16", "LT-18", "LT-19", "LT-2",
           "LT-20", "LT-21", "LT-24", "LT-25", "LT-26",
           "LT-27", "LT-28", "LT-29", "LT-30", "LT-4",
           "LT-5", "LT-7", "LT-9", "LT-23")

lt<-as.factor(c(male,female))

##extract year

fish$YEAR<-year(fish$DATETIME)
fish$MONTH<-month(fish$DATETIME)
fish$DAY<-day(fish$DATETIME)
fish$HOUR<-hour(fish$DATETIME)

#assign sex M = male, F = female, U = unknown
fish$SEX<-"U"
fish$SEX[fish$TRANSMITTER %in% male]<-"M"
fish$SEX[fish$TRANSMITTER %in% female]<-"F"



#shore distance
require(SDMTools)
shore.dist<-read.asc("~/Documents/Masters Thesis/Data/GIS/Alexie Lake/Alexie Raster/shoredist.asc")
#assign raster class object
shore.dist<-raster.from.asc(shore.dist)
#assign projection
projection(shore.dist)<- "+proj=utm +zone=11 +ellps=WGS84"
#extract distance from shore for each point
fish$SHOREDIST<-extract(shore.dist, fish[,c(1,2)])
fish$SHOREDIST[is.na(fish$SHOREDIST)]<-0

fish$BOTTOMDEPTH<-(extract(r,fish[c("UTM.X","UTM.Y")]))
fish$OFFBOTTOM<-fish$BOTTOMDEPTH-fish$DEPTH

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

fish$SEASON<-"SUMMER"
fish$SEASON[fish$DATETIME>=as.POSIXct("2013-08-30",tz="MST") & 
                 fish$DATETIME<as.POSIXct("2013-10-01",tz="MST")]<-"SPAWN"
fish$SEASON[fish$DATETIME>=as.POSIXct("2013-10-01",tz="MST")]<-"POSTSPAWN"
fish$SEASON[fish$DATETIME>=as.POSIXct("2013-05-15",tz="MST") & 
                   fish$DATETIME<as.POSIXct("2013-06-25",tz="MST")]<-"SPRING"


clust.param<-data.frame(Season=NULL,
                        ID=NULL,
                        MinPts=NULL,
                        EpsSpace=NULL,
                        EpsTime=NULL)
season<-unique(fish$SEASON)
tag<-c("LT-31","LT-32","LT-33","LT-34","LT-35","LT-37","LT-38","LT-40","LT-41","LT-43","LT-44")

for(i in 1:length(season)){
  
  fishy<-fish[fish$SEASON==season[i],]
  #tag<-unique(fishy$TRANSMITTER)
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
    
    clust.p<-data.frame(Season=season[i],
                        ID=tag[ii],
                        MinPts=k,
                        EpsSpace=hd,
                        EpsTime=ht)
    
    clust.param<-rbind(clust.param,clust.p)
  }
}
rm(fishy)
clust.param

                
#for each parameter a take the rounded mean of from the sampled animals

#Spring
min_pts_spring<-round(summary(clust.param$MinPts[clust.param$Season=="SPRING"]))
eps_space_spring<-round(summary(clust.param$EpsSpace[clust.param$Season=="SPRING"]))
eps_time_spring<-round(summary(clust.param$EpsTime[clust.param$Season=="SPRING"]))

#Summer
min_pts_summer<-round(summary(clust.param$MinPts[clust.param$Season=="SUMMER"]))
eps_space_summer<-round(summary(clust.param$EpsSpace[clust.param$Season=="SUMMER"]))
eps_time_summer<-round(summary(clust.param$EpsTime[clust.param$Season=="SUMMER"]))


#Spawn
min_pts_spawn<-round(summary(clust.param$MinPts[clust.param$Season=="SPAWN"]))
eps_space_spawn<-round(summary(clust.param$EpsSpace[clust.param$Season=="SPAWN"]))
eps_time_spawn<-round(summary(clust.param$EpsTime[clust.param$Season=="SPAWN"]))


cluster<-data.frame(Season=season,
                    Min_Pts_Mean=c(min_pts_spring[["Mean"]],min_pts_summer[["Mean"]],min_pts_spawn[["Mean"]]),
                    Min_Pts_Min=c(min_pts_spring[["Min."]],min_pts_summer[["Min."]],min_pts_spawn[["Min."]]),
                    Min_Pts_Max=c(min_pts_spring[["Max."]],min_pts_summer[["Max."]],min_pts_spawn[["Max."]]),
                    
                    Eps_Space_Mean=c(eps_space_spring[["Mean"]],eps_space_summer[["Mean"]],eps_space_spawn[["Mean"]]),
                    Eps_Space_Min=c(eps_space_spring[["Min."]],eps_space_summer[["Min."]],eps_space_spawn[["Min."]]),
                    Eps_Space_Max=c(eps_space_spring[["Max."]],eps_space_summer[["Max."]],eps_space_spawn[["Max."]]),
                    
                    Eps_Time_Mean=c(eps_time_spring[["Mean"]],eps_time_summer[["Mean"]],eps_time_spawn[["Mean"]]),
                    Eps_Time_Min=c(eps_time_spring[["Min."]],eps_time_summer[["Min."]],eps_time_spawn[["Min."]]),
                    Eps_Time_Max=c(eps_time_spring[["Max."]],eps_time_summer[["Max."]],eps_time_spawn[["Max."]])
)

min_pts<-min_pts_spawn[["Mean"]]
eps_space<-eps_space_spawn[["Mean"]]
eps_time<-eps_time_spawn[["Mean"]]


for(ii in 1:length(tag)){
  fishi<-fish[fish$TRANSMITTER==tag[ii],]
  fishi<-fishi[order(fishi$DATETIME),]
  x=fishi$UTM.X
  y=fishi$UTM.Y
  t=fishi$DATETIME
  clust<-stdbscan(lon=x, lat=y, time=t, eps_space=eps_space, eps_time=eps_time, deltaE=2,SD=TRUE,MinPts=min_pts)
  
  print(clust)
  print(paste0("Progress = ", round((ii/length(tag))*100),"%"))
  print(tag[ii])
  fishi$CLUSTER<-clust$cluster
  assign(paste("LT",substr(tag[ii],start=4,stop=5),sep=""),fishi)
  assign(paste("clusterLT",substr(tag[ii],start=4,stop=5),sep=""),clust)
}

fishlist<-list(LT31,LT32,LT33,LT34,LT35,LT37,LT38,LT40,LT41,LT43,LT44)
rm(LT31,LT32,LT33,LT34,LT35,LT37,LT38,LT40,LT41,LT43,LT44)
rm(fishi)

unique(fish$TRANSMITTER)
#need spawn_pts and egg_pts from All Spawn Sites.R

#spawnpoly<-gBuffer(spawn_pts,byid=TRUE,width=20)

for(i in 1:length(fishlist)){
  spawnclust=NULL
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))
  plot(shore_outline)
  title(tag[i])
  
  lines(fishi$UTM.X,fishi$UTM.Y,col=rgb(0,0,0,0.2),pch=19,cex=0.2)
  points(xy,col=rgb(0,0,0,0.2),pch=19,cex=0.02)
  
  for(ii in 1:max(fishi$CLUSTER)){
    new.xy<-xy[fishi$CLUSTER==ii]
    #clusters must have a median point distance < 40 m from shore and median depth < 5 m
    #if(median(fishi$SHOREDIST[fishi$CLUSTER==ii],na.rm=TRUE)<40 & median(fishi$DEPTH[fishi$CLUSTER==ii],na.rm=TRUE)<5){
    #if(length(new.xy)>=8){
    colour<-"red"
    cp=mcp(new.xy,percent=100)
    #cp1=gBuffer(cp, width=50.0, quadsegs=5, capStyle="SQUARE",
    #       joinStyle="MITRE", mitreLimit=1.0)
    #plot(cp1,add=T,col=colour+1)
    plot(cp,add=TRUE,col=colour)
    spawnclust2<-data.frame(CLUSTER=ii,
                            X=mean(cp@bbox[1,c(1,2)]),
                            Y=mean(cp@bbox[2,c(1,2)]),
                            AREA=gArea(cp),
                            SPAWNSITE=ifelse(any(gOverlaps(cp,spawnpoly, byid=TRUE)),1,0),
                            EGGS=ifelse(any(gOverlaps(cp,spawnpoly[spawnpoly$EGGS==1,], byid=TRUE)),1,0))
    spawnclust<-rbind(spawnclust,spawnclust2)
    
    #}
    #}
  }
  assign(paste("spawnclustLT",substr(tag[i],start=4,stop=5),sep=""),spawnclust)
}


#cluster list; NOTE:LT32 never spawns
clustlist<-list(spawnclustLT31,spawnclustLT32,spawnclustLT33,spawnclustLT34,spawnclustLT35,spawnclustLT37,spawnclustLT38,spawnclustLT40,spawnclustLT41,spawnclustLT43,spawnclustLT44)
rm(spawnclustLT31,spawnclustLT32,spawnclustLT33,spawnclustLT34,spawnclustLT35,spawnclustLT37,spawnclustLT38,spawnclustLT40,spawnclustLT41,spawnclustLT43,spawnclustLT44)

##### param graphs ####
#split fish into months
#select spring data 2 days after the last surgery
fishspring<-fish[fish$SEASON=="SPRING" & fish$DAY>15,]
#select summer data
fishsummer<-fish[fish$SEASON=="SUMMER",]

#select spawn data (no october)
fishspawn<-fish[fish$SEASON=="SPAWN" & fish$MONTH<10,]

fishspawn$SEX<-as.factor(fishspawn$SEX)

save(fish,fishspawn,shore_outline,depth15,file="fishspawn.RData")

#Spawn
quartz()
plot(shore_outline,main="Lake trout spawn, 2013")
points(UTM.Y~UTM.X,data=fishspawn[fishspawn$SEX=="F",],col=rgb(1,0,0,0.1),pch=19,cex=0.1)
points(UTM.Y~UTM.X,data=fishspawn[fishspawn$SEX=="M",],col=rgb(0,0,1,0.1),pch=19,cex=0.1)

acc<-fishspawn$ACCEL[!is.na(fishspawn$ACCEL)]
sacc<-ltspawnclust3$ACCEL
ltspawnclust3$SEX<-as.factor(ltspawnclust3$SEX)
ltspawnclust3$SPAWN<-as.factor(ltspawnclust3$SPAWN)


shapiro.test(ltspawnclust3$ACCEL)
bartlett.test(ltspawnclust3$ACCEL,ltspawnclust3$SPAWN)
t.test(ltspawnclust3$ACCEL~ltspawnclust3$SEX)

sex<-unique(fishspawn$SEX)
means<-c(mean(fishspawn$ACCEL[!is.na(fishspawn$ACCEL) & fishspawn$SEX=="M"]),mean(fishspawn$ACCEL[!is.na(fishspawn$ACCEL) & fishspawn$SEX=="F"]))
sds<-c(sd(fishspawn$ACCEL[!is.na(fishspawn$ACCEL) & fishspawn$SEX=="M"]),sd(fishspawn$ACCEL[!is.na(fishspawn$ACCEL) & fishspawn$SEX=="F"]))
acc<-data.frame(SEX=as.factor(sex),
                MEAN=means,
                SD=sds
)


plot(acc$SEX,acc$MEAN,type="h")

arrows(MEAN+SD,mean,data=fishspawn[!is.na(fishspawn$ACCEL),])


chisq.test(table(ltspawnclust3$SEX))

hist(log(acc))

ks.test(log(acc),"pnorm")
bartlett.test(log(acc),fishspawn$SEX[!is.na(fishspawn$ACCEL)])

#t.test(acc~fishspawn$SEX[!is.na(fishspawn$ACCEL)])
wilcox.test(acc~fishspawn$SEX[!is.na(fishspawn$ACCEL)])

#Nearshore


## calculate the density - don't plot yet

densALL<-density(fishspawn$SHOREDIST, from=0)
densFemale_Spawn <- density(fishspawn$SHOREDIST[fishspawn$SEX=="F"],from=0, bw = densALL$bw)
densMale_Spawn<- density(fishspawn$SHOREDIST[fishspawn$SEX=="M"],from=0, bw = densALL$bw)


## calculate the range of the graph
xlim <- range(0,450)
ylim <- range(0,0.012)
#pick the colours
FemaleCol <- rgb(1,0,0,0.5)
MaleCol <- rgb(0,0,1,0.5)
## plot the carrots and set up most of the plot parameters
plot(densFemale_Spawn, xlim = xlim, ylim = ylim, 
     xlab = 'Distance to Shore',
     main="",
     type="n",
     bty="n",
     axes=FALSE
)

#put our density plots in
polygon(c(densMale_Spawn$x,rev(densMale_Spawn$x)),c(densMale_Spawn$y,rep(0,512)), density = -1, col = MaleCol, border = "blue")
polygon(c(densFemale_Spawn$x,rev(densFemale_Spawn$x)),c(densFemale_Spawn$y,rep(0,512)), density = -1, col = FemaleCol, border = "red")



# Draw our own x-axis
axis(1,pos=0)
axis(1,at=xlim,labels=c("",""),pos=0,lwd.ticks=0)

axis(2,las=1,pos=0)


## add a legend in the corner
legend("topright",c('Male', "Female"),
       fill = c(MaleCol,FemaleCol), bty = 'n',
       border = NA)

abline(v=50)
mtext("50",side=1,line=0.5,at=50)

#Depth 

quartz()
par(mfrow=c(2,1))
par(mar=c(0,4.1,4.1,2.1))
#5.1 4.1 4.1 2.1
#males
plot(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15",],
     main="",
     xlab="",
     ylab="",
     type="n",
     bty="n",
     axes=FALSE,
     ylim=c(32,0))

#axis.POSIXct(1,at=c(as.POSIXct("2013-08-15", tz="MST"),as.POSIXct("2013-09-01",tz="MST") ,as.POSIXct("2013-09-15", tz="MST"),as.POSIXct("2013-10-01", tz="MST"),as.POSIXct("2013-10-15", tz="MST")), format="%d %b")

axis(2,
     las=1
)
     #pos=as.numeric(as.POSIXct("2013-08-15", tz="MST")))
box()




points(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15" & fish$SEX=="M",],
       col=rgb(0,0,1,0.1),pch=15)
lines(depth15$Date,depth15$Depth,lwd=4)
text(as.numeric(as.POSIXct("2013-08-14", tz="MST")),0.5,"(A)", font=2)

#females
par(mar=c(5.1,4.1,0,2.1))
plot(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15",],
     main="",
     xlab="",
     ylab="",
     type="n",
     bty="n",
     axes=FALSE,
     ylim=c(32,0))

axis.POSIXct(1,at=c(as.POSIXct("2013-08-15", tz="MST"),as.POSIXct("2013-09-01",tz="MST") ,as.POSIXct("2013-09-15", tz="MST"),as.POSIXct("2013-10-01", tz="MST"),as.POSIXct("2013-10-15", tz="MST")),
             format="%d %b")

axis(2,
     las=1)
box()

points(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15" & fish$SEX=="F",],
       col=rgb(1,0,0,0.1),pch=15)
lines(depth15$Date,depth15$Depth,lwd=4)
text(as.numeric(as.POSIXct("2013-08-14", tz="MST")),0.5,"(B)", font=2)
mtext("Depth (m)", side=2, line=3, adj=1.2)


#Actual Depth

library(ks)

## calculate the density - don't plot yet
densFemale_Spawn <- density(fishspawn$DEPTH[fishspawn$SEX=="F"],from=0, na.rm=TRUE)
densMale_Spawn<- density(fishspawn$DEPTH[fishspawn$SEX=="M"],from=0, na.rm=TRUE)

densFemale_Spawn$bw
densFemale_Spawn$n

densMale_Spawn$bw
densMale_Spawn$n

## calculate the range of the graph
xlim <- range(0,32)
ylim <- range(0,0.2)
#pick the colours
FemaleCol <- rgb(1,0,0,0.5)
MaleCol <- rgb(0,0,1,0.5)
## plot the carrots and set up most of the plot parameters
plot(densMale_Spawn, xlim = xlim, ylim = ylim, 
     xlab = 'Depth',
     main="",
     type="n",
     bty="n",
     axes=FALSE
)

#put our density plots in
polygon(c(densMale_Spawn$x,rev(densMale_Spawn$x)),c(densMale_Spawn$y,rep(0,512)), density = -1, col = MaleCol, border = "blue")
polygon(c(densFemale_Spawn$x,rev(densFemale_Spawn$x)),c(densFemale_Spawn$y,rep(0,512)), density = -1, col = FemaleCol, border = "red")


# Draw our own x-axis
axis(1,pos=0)
axis(1,at=xlim,labels=c("",""),pos=0,lwd.ticks=0)

axis(2,las=1,pos=0)


## add a legend in the corner
legend("topright",c('Male', "Female"),
       fill = c(MaleCol,FemaleCol), bty = 'n',
       border = NA)

abline(v=4)
mtext("4",side=1,line=0.5,at=4)

require(mixtools)
fit<-normalmixEM(fishspawn$DEPTH[fishspawn$SEX=="M" & !is.na(fishspawn$DEPTH)],lambda=0.5,mu=c(2,12),sigma=5)

#Depth Contour

## calculate the density - don't plot yet
densFemale_Spawn <- density(fishspawn$BOTTOMDEPTH[fishspawn$SEX=="F"],from=0)
densMale_Spawn<- density(fishspawn$BOTTOMDEPTH[fishspawn$SEX=="M"],from=0)

## calculate the range of the graph
xlim <- range(0,32)
ylim <- range(0,0.1)
#pick the colours
FemaleCol <- rgb(1,0,0,0.5)
MaleCol <- rgb(0,0,1,0.5)
## plot the carrots and set up most of the plot parameters
plot(densMale_Spawn, xlim = xlim, ylim = ylim, 
     xlab = 'Depth Contour',
     main = "Lake Trout Spawn, 2013",
     type="n",
     bty="n",
     axes=FALSE
)

#put our density plots in
polygon(c(densMale_Spawn$x,rev(densMale_Spawn$x)),c(densMale_Spawn$y,rep(0,512)), density = -1, col = MaleCol, border = "blue")
polygon(c(densFemale_Spawn$x,rev(densFemale_Spawn$x)),c(densFemale_Spawn$y,rep(0,512)), density = -1, col = FemaleCol, border = "red")


# Draw our own x-axis
axis(1,pos=0)
axis(1,at=xlim,labels=c("",""),pos=0,lwd.ticks=0)

axis(2,las=1,pos=0)


## add a legend in the corner
legend("topright",c('Male', "Female"),
       fill = c(MaleCol,FemaleCol), bty = 'n',
       border = NA)

abline(v=4)
mtext("4",side=1,line=0.5,at=4)



#Depth Above Bottom

## calculate the density - don't plot yet
densFemale_Spawn<- density(fishspawn$OFFBOTTOM[fishspawn$SEX=="F"],from=0,to=32,na.rm=TRUE)

densMale_Spawn<-density(fishspawn$OFFBOTTOM[fishspawn$SEX=="M"],from=0,to=32,na.rm=TRUE)
## calculate

## calculate the range of the graph
xlim <- range(0,32)
ylim <- range(0,0.15)
#pick the colours
FemaleCol <- rgb(1,0,0,0.5)
MaleCol <- rgb(0,0,1,0.5)
## plot the carrots and set up most of the plot parameters
plot(densMale_Spawn, xlim = xlim, ylim = ylim, 
     xlab = 'Depth Above Bottom',
     main = "Lake Trout Spawn, 2013",
     type="n",
     bty="n",
     axes=FALSE
)

#put our density plots in
polygon(c(densMale_Spawn$x,rev(densMale_Spawn$x)),c(densMale_Spawn$y,rep(0,512)), density = -1, col = MaleCol, border = "blue")
polygon(c(densFemale_Spawn$x,rev(densFemale_Spawn$x)),c(densFemale_Spawn$y,rep(0,512)), density = -1, col = FemaleCol, border = "red")


# Draw our own x-axis
axis(1,pos=0)
axis(1,at=xlim,labels=c("",""),pos=0,lwd.ticks=0)

axis(2,las=1,pos=0)


## add a legend in the corner
legend("topright",c('Male', "Female"),
       fill = c(MaleCol,FemaleCol), bty = 'n',
       border = NA)

abline(v=3)
mtext("3",side=1,line=0.5,at=3)

#####



ltspawnclust<-NULL

for(i in 1:length(fishlist)){
  fishi=fishlist[[i]]
  
  spawnclust=clustlist[[i]]
  if(length(spawnclust$CLUSTER)>0){
    fishi<-fishi[fishi$CLUSTER %in% spawnclust$CLUSTER,]
    
    
    for(ii in 1: length(spawnclust$CLUSTER)){
      fishi2<-fishi[fishi$CLUSTER==spawnclust$CLUSTER[ii],]
      ltspawnclust2<-data.frame(TRANSMITTER=fishi2$TRANSMITTER[1],
                                CLUSTER=spawnclust$CLUSTER[ii],
                                CLUST.X=spawnclust$X[ii],
                                CLUST.Y=spawnclust$Y[ii],
                                AREA=spawnclust$AREA[ii],
                                SPAWNSITE=spawnclust$SPAWNSITE[ii],
                                EGGS=spawnclust$EGGS[ii],
                                STARTTIME=min(fishi2$DATETIME),
                                TOTALTIME=as.numeric(difftime(max(fishi2$DATETIME),min(fishi2$DATETIME),units="mins")),
                                DEPTH=mean(fishi2$DEPTH,na.rm=TRUE),
                                MINDEPTH=min(fishi2$DEPTH,na.rm=TRUE),
                                MAXDEPTH=max(fishi2$DEPTH,na.rm=TRUE),
                                BOTTOMDEPTH=mean(fishi2$BOTTOMDEPTH,na.rm=TRUE),
                                PERCUNDER4M=length(fishi2$BOTTOMDEPTH[fishi2$BOTTOMDEPTH<4])/length(fishi2$BOTTOMDEPTH),
                                MAXBOTTOMDEPTH=max(fishi2$BOTTOMDEPTH,na.rm=TRUE),
                                MINBOTTOMDEPTH=min(fishi2$BOTTOMDEPTH,na.rm=TRUE),
                                SHOREDIST=mean(fishi2$SHOREDIST,na.rm=TRUE),
                                PERCUNDER50M=length(fishi2$SHOREDIST[fishi2$SHOREDIST<50])/length(fishi2$SHOREDIST),
                                ACCEL=mean(fishi2$ACCEL,na.rm=TRUE),
                                MAXACCEL=max(fishi2$ACCEL,na.rm=TRUE),
                                PERCOVERACCELCUT=length(fishi2$ACCEL[fishi2$ACCEL>1.3 & !is.na(fishi2$ACCEL)])/length(fishi2$ACCEL[!is.na(fishi2$ACCEL)])
      )
      ltspawnclust<-rbind(ltspawnclust,ltspawnclust2)
    }
    assign(paste("spawnclustLT",substr(tag[i],start=4,stop=5),sep=""),fishi)
  }
}

ltspawnclust$SEX<-"U"
ltspawnclust$SEX[ltspawnclust$TRANSMITTER %in% male]<-"M"
ltspawnclust$SEX[ltspawnclust$TRANSMITTER %in% female]<-"F"

plot(shore_outline,main="All Clusters")
points(CLUST.Y~CLUST.X,data=ltspawnclust,pch=19,cex=0.3)
#area_cut=mean(ltspawnclust$AREA)+2*sd(ltspawnclust$AREA)
#data1=ltspawnclust[ltspawnclust$AREA<area_cut,]


plot(shore_outline,main="Potenital Spawn Clusters")
data1=ltspawnclust[ltspawnclust$PERCUNDER50M>.50,]
points(CLUST.Y~CLUST.X,data=data1,pch=19,col=rgb(1,0,0,0.3))
data2=data1[data1$PERCUNDER4M>.50,]
points(CLUST.Y~CLUST.X,data=data2,pch=19,cex=0.5,col=rgb(0,0,1,1))


ltspawnclust2<-data2
ltspawnclust2$SITE<-NA

#same site spawning
ltspawnclust3=NULL
ctags<-unique(ltspawnclust2$TRANSMITTER)
for(i in 1:length(ctags)){
  clust<-ltspawnclust2[ltspawnclust2$TRANSMITTER==ctags[i],]
  clusters<-unique(clust$CLUSTER)
  for(ii in 1:length(clusters)){
    clust2<-clust[clust$CLUSTER==clusters[ii],]
    if(is.na(clust2$SITE)){
      d<-pointDistance(clust2[,c("CLUST.X","CLUST.Y")],clust[,c("CLUST.X","CLUST.Y")],allpairs=TRUE,longlat=FALSE)
      #the centre of clusters must be 70 m apart to be concidered "unique sites"
      clust$SITE[d<70]<-ii
    }
  }
  ltspawnclust3=rbind(clust,ltspawnclust3)
}

ltspawnclust3$HOUR<-hour(ltspawnclust3$STARTTIME)+minute(ltspawnclust3$STARTTIME)/60+second(ltspawnclust3$STARTTIME)/3600
ltspawnclust3$DAY<-day(ltspawnclust3$STARTTIME)
rownames(ltspawnclust3)<-seq(1,dim(ltspawnclust3)[1],by=1)

fishbind<-do.call(rbind,fishlist)

#mean ACCEL+2 sd =accel of 1.3
accel_cut<-mean(fishbind$ACCEL[fishbind$DATETIME>="2013-08-30" &fishbind$DATETIME<"2013-10-01" & fishbind$CLUSTER!=0],na.rm=TRUE)+2*(sd(fishbind$ACCEL[fishbind$DATETIME>="2013-08-30" &fishbind$DATETIME<"2013-10-01" & fishbind$CLUSTER!=0],na.rm=TRUE))
ltspawnclust3$SPAWN<-0
ltspawnclust3$SPAWN[ltspawnclust3$MAXACCEL>accel_cut]<-1

#plots 6.15 x 4.76
plot(shore_outline)
points(CLUST.Y~CLUST.X,data=ltspawnclust3)
points(CLUST.Y~CLUST.X,data=ltspawnclust3[ltspawnclust3$SPAWN==1,],pch=19,col=rgb(1,0,0,1))
points(CLUST.Y~CLUST.X,data=ltspawnclust3[ltspawnclust3$SPAWNSITE==1,],pch=19,col=rgb(0,1,0,1))
points(CLUST.Y~CLUST.X,data=ltspawnclust3[ltspawnclust3$EGGS==1,],pch=19,col=rgb(0,0,1,1))



points(UTM.Y~UTM.X,data=spawn_pts,pch=19,cex=0.5,col="blue")
points(UTM.Y~UTM.X,data=spawn_pts[spawn_pts$EGGS==1,],pch=19,cex=0.5,col="green")

ltspawnclust3$DATE<-format(ltspawnclust3$STARTTIME,"%Y-%m-%d")

#all clusters
spawn.num=NULL
spawn.num$DATE<-as.POSIXct(names(table(ltspawnclust3$DATE)),format="%Y-%m-%d", tz="MST")
spawn.num$FREQ<-as.numeric(table(ltspawnclust3$DATE))
#all spawn clusters
spawn.num2=NULL
spawn.num2$SPAWN<-as.numeric(table(ltspawnclust3$DATE[ltspawnclust3$SPAWN==1]))
spawn.num2$DATE<-as.POSIXct(names(table(ltspawnclust3$DATE[ltspawnclust3$SPAWN==1])),format="%Y-%m-%d", tz="MST")
spawn.num<-as.data.frame(spawn.num)
spawn.num2<-as.data.frame(spawn.num2)
spawn.num<-merge(spawn.num,spawn.num2,by="DATE",all=TRUE)

#male freq
spawn.num2=NULL
spawn.num2$DATE<-as.POSIXct(names(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="M"])),format="%Y-%m-%d", tz="MST")
spawn.num2$FREQ_MALE<-as.numeric(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="M"]))
spawn.num2<-as.data.frame(spawn.num2)
spawn.num<-merge(spawn.num,spawn.num2,by="DATE",all=TRUE)
#female freq
spawn.num2=NULL
spawn.num2$DATE<-as.POSIXct(names(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="F"])),format="%Y-%m-%d", tz="MST")
spawn.num2$FREQ_FEMALE<-as.numeric(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="F"]))
spawn.num2<-as.data.frame(spawn.num2)
spawn.num<-merge(spawn.num,spawn.num2,by="DATE",all=TRUE)

#male spawn
spawn.num2=NULL
spawn.num2$DATE<-as.POSIXct(names(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="M" & ltspawnclust3$SPAWN==1])),format="%Y-%m-%d", tz="MST")
spawn.num2$SPAWN_MALE<-as.numeric(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="M" & ltspawnclust3$SPAWN==1]))
spawn.num2<-as.data.frame(spawn.num2)
spawn.num<-merge(spawn.num,spawn.num2,by="DATE",all=TRUE)
#female spawn
spawn.num2=NULL
spawn.num2$DATE<-as.POSIXct(names(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="F" & ltspawnclust3$SPAWN==1])),format="%Y-%m-%d", tz="MST")
spawn.num2$SPAWN_FEMALE<-as.numeric(table(ltspawnclust3$DATE[ltspawnclust3$SEX=="F" & ltspawnclust3$SPAWN==1]))
spawn.num2<-as.data.frame(spawn.num2)
spawn.num<-merge(spawn.num,spawn.num2,by="DATE",all=TRUE)

spawn.num[is.na(spawn.num)]<-0



dateseq= seq.POSIXt(from=as.POSIXct("2013-08-15",tz="MST"),to=as.POSIXct("2013-11-01", tz="MST"),by="day")


save(ltspawnclust3,temp,spawn.num,dateseq,suntimes,fishlist,file="SpawnClust.RData")

load("SpawnClust.RData")
quartz()
#par(mfrow=c(2,1))
#par(mar=c(0,2.1,4.1,4.1))

peakstart=min(spawn.num$DATE[spawn.num$FREQ>=(max(spawn.num$FREQ)/2)])
peakend=max(spawn.num$DATE[spawn.num$FREQ>=(max(spawn.num$FREQ)/2)])
spawnstart=min(spawn.num$DATE)
spawnend=max(spawn.num$DATE)


plot(FREQ~DATE,data=spawn.num,col=rgb(0,0,0,0.3),ylim=c(0,25),type="h",xlim=c(min(dateseq),max(dateseq)),main="",ylab="",xlab="",axes=FALSE,lwd=7,lend=1)

lines(SPAWN~DATE,data=spawn.num,col=rgb(0,0,0,0.7),type="h",lwd=7,lend=1)

lines(c(peakstart-12*60*60,peakstart-12*60*60),c(0,25))
lines(c(peakend+12*60*60,peakend+12*60*60),c(0,25))

axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d",pos=0)
axis(2,pos=min(dateseq))


lines(Temp~Date,data=temp[temp$Depth==0 & temp$Date >="2013-08-15" & temp$Date <="2013-11-01",],lwd=2)


#quartz()
#par(mfrow=c(2,1))
#par(mar=c(0,4.1,4.1,2.1))

plot(FREQ~DATE,data=spawn.num,col=rgb(0,0,0,0.3),ylim=c(-15,15),type="n",xlim=c(min(dateseq),max(dateseq)),main="",ylab="Frequency",xlab="",axes=FALSE)

lines(FREQ_MALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(0,0,1,0.3))
lines(SPAWN_MALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(0,0,1,0.7))
lines(-FREQ_FEMALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(1,0,0,0.3))
lines(-SPAWN_FEMALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(1,0,0,0.7))

#lines(c(peakstart-12*60*60,peakstart-12*60*60),c(-10,15))
#lines(c(peakend+12*60*60,peakend+12*60*60),c(-10,15))

axis.POSIXct(1,at=c(dateseq[1],dateseq[79]),labels=FALSE,lwd.ticks=0,pos=0)

ylabs <- pretty(c(-10,15))
axis(2, at = ylabs, labels = abs(ylabs), pos=min(dateseq),las=1)


axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d",pos=-10)

unique(ltspawnclust3$TRANSMITTER)


#axis.POSIXct(1,at=seq(1,62,length.out=5),seq.POSIXt(from=min(dateseq),to=max(dateseq),length.out=5),format="%b-%d",pos=0)

#axis(2,las=1,pos=min(dateseq))


#####suntimes
require(maptools)
#create spatial point
point<-fish[1,c("LAT","LON")]
coordinates(point)<- ~LON+LAT
proj4string(point)<-CRS("+proj=longlat +datum=WGS84")

#date sequence
#dateseq= seq.POSIXt(from=as.POSIXct("2013-08-15",tz="MST"),to=as.POSIXct("2013-11-01", tz="MST"),by="day")

#sun times data frame
suntimes<-data.frame(DATE=dateseq,
                     SUNRISE=crepuscule(point, dateseq, solarDep = 0, direction="dawn",POSIXct.out=TRUE)$time,
                     DAWN=crepuscule(point, dateseq, solarDep = 6, direction="dawn",
                                     POSIXct.out=TRUE)$time,
                     DUSK=crepuscule(point, dateseq, solarDep = 6, direction= "dusk",
                                     POSIXct.out=TRUE)$time,
                     SUNSET=crepuscule(point, dateseq, solarDep = 0, direction="dusk",POSIXct.out=TRUE)$time
)

suntimes$SUNRISE<-hour(suntimes$SUNRISE)+minute(suntimes$SUNRISE)/60+second(suntimes$SUNRISE)/3600

suntimes$SUNSET<-hour(suntimes$SUNSET)+minute(suntimes$SUNSET)/60+second(suntimes$SUNSET)/3600

suntimes$DAWN<-hour(suntimes$DAWN)+minute(suntimes$DAWN)/60+second(suntimes$DAWN)/3600

suntimes$DUSK<-hour(suntimes$DUSK)+minute(suntimes$DUSK)/60+second(suntimes$DUSK)/3600

ltspawnclust3$DATE<-as.POSIXct(ltspawnclust3$DATE,tz="MST")
ltspawnclust3$TIME<-NA
for(i in 1:dim(ltspawnclust3)[1]){
  times<-suntimes[suntimes$DATE==ltspawnclust3$DATE[i],]
  ltspawnclust3$TIME[i]<-ifelse(ltspawnclust3$HOUR[i] > times$SUNSET | ltspawnclust3$HOUR[i] < times$SUNRISE, "NIGHT","DAY")
}
  


ltspawnclust3$HOURPLOT<-hour(ltspawnclust3$STARTTIME)
ltspawnclust3$HOURPLOT[ltspawnclust3$HOURPLOT==0]<-ltspawnclust3$HOURPLOT[ltspawnclust3$HOURPLOT==0]+24

require(psych)
circadian.mean(ltspawnclust3$HOUR)
circadian.sd(ltspawnclust3$HOUR)

y <- table(factor(ltspawnclust3$HOURPLOT, levels=1:24))
y<-as.data.frame(y)
y$Var1<-as.numeric(y$Var1)-1

y.smooth<-loess(Freq~Var1,data=y)

y.smooth<-predict(y.smooth, se=T)

y <- table(factor(ltspawnclust3$HOURPLOT[ltspawnclust3$SPAWN==1], levels=1:24))
y<-as.data.frame(y)
y$Var1<-as.numeric(y$Var1)-1

y.smooth2<-loess(Freq~Var1,data=y)

y.smooth2<-predict(y.smooth2, se=T)


quartz()
layout(matrix(c(1,1,1,2,1,1,1,2),2,4,byrow=TRUE))
#par(mfrow=c(1,2))
par(mar=c(5.1,4.1,4.1,0),
    oma=c(2,2,0,0))

plot(HOURPLOT~STARTTIME, data=ltspawnclust3, col=rgb(0,0,0,0),pch=19, xlab="Date",ylab="Hour", xlim=c(min(dateseq),max(dateseq)),ylim=c(0.5,24.5), axes=FALSE,xaxs="i",yaxs="i")



axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d",pos=0.5)


axis(2,at=seq(0,24,by=4),
     las=1)

#sunrise
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$SUNRISE,rep(0.5,dim(suntimes)[1])), density = -1, col = "grey",border="grey")
#dawn
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$DAWN,rep(0.5,dim(suntimes)[1])), density = -1, col = "darkgrey", border = "darkgrey")


#sunset
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$SUNSET,rep(24.5,dim(suntimes)[1])), density = -1, col = "grey",border="grey")
#dusk
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$DUSK,rep(24.5,dim(suntimes)[1])), density = -1, col = "darkgrey", border = "darkgrey")

points(HOURPLOT~STARTTIME, data=ltspawnclust3,col=rgb(1,0,0,1))
#points(HOUR~STARTTIME, data=ltspawnclust3[ltspawnclust3$SPAWN==0,], pch=19, col=rgb(0,0,1,1))
points(HOURPLOT~STARTTIME, data=ltspawnclust3[ltspawnclust3$SPAWN==1,],pch=19,col=rgb(1,0,0,1))


box(lwd=2)


par(mar=c(5.1,0,4.1,1))


plot(y.smooth$fit,y$Var1,pch=19,ylim=c(0.5,24.5), xlim=c(0,16),xlab="", ylab="",type="n",axes=FALSE,xaxs="i",yaxs="i")
axis(1,at=c(0,16),labels=FALSE,lwd.ticks=0,pos=0.5,lwd=2)
axis(2,at=c(0.5,24.5),labels=FALSE,lwd.ticks=0,pos=0,lwd=2)
axis(3,at=seq(0,16,4), pos=24.5,lwd=2)
axis(4,at=c(0.5,24.5),labels=FALSE,lwd.ticks=0,pos=16,lwd=2)
mtext("Count",side=3,line=3,cex=0.6)
#box()


#polygon(c(y.smooth$fit+ qt(0.975,y.smooth$df)*y.smooth$se,rev(y.smooth$fit- qt(0.975,y.smooth$df)*y.smooth$se)), c(y$Var1,rev(y$Var1)), density = -1, col = rgb(0,0,0,0.2),border=FALSE)

arrows(x0=y.smooth$fit+ qt(0.975,y.smooth$df)*y.smooth$se,
       y0=y$Var1,
       x1=y.smooth$fit- qt(0.975,y.smooth$df)*y.smooth$se, 
       angle=90,code=3,length=0.05)


arrows(x0=y.smooth2$fit+ qt(0.975,y.smooth2$df)*y.smooth2$se,
       y0=y$Var1,
       x1=y.smooth2$fit- qt(0.975,y.smooth2$df)*y.smooth2$se, 
       angle=90,code=3,length=0.05)

points( y.smooth$fit,y$Var1,pch=19, col="white")
lines( y.smooth$fit,y$Var1,type="b",lwd=2, col="red")
lines( y.smooth2$fit,y$Var1,type="b",lwd=2, pch=19,col="red")

#points( y.smooth$y[x.mid+offset]^2, x.mid-1, col="red")
#points( y.smooth2$y[x.mid+offset]^2, x.mid-1,pch=19, col="red")
#sunrise
#polygon(c(y.smooth$y[x.mid+offset]^2,rep(0,length(y.smooth$y[x.mid+offset]^2))),c(x.mid-1,rev(x.mid-1)), density = -1, col = rgb(1,0,0,0.5),border=FALSE)


plot(spawn$DATE,spawn$TRANSMITTER,pch=15, cex=0.7,axes=FALSE,xlim=c(min(dateseq),max(dateseq)),xlab="Fish ID",ylab="")
axis(2,at=c(1:length(tag)),labels=tag,las=1)
axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d")
box(lwd=2)

##daily displacement


fishbind$DATE<-format(fishbind$DATETIME,"%Y-%m-%d")
for(i in 1:length(tag)){
data.unclass<-fishbind[fishbind$TRANSMITTER==tag[i] & fishbind$SEASON=="SPAWN",]
x<-data.unclass$UTM.X
y<-data.unclass$UTM.Y
#add x and y data as real and imaginary numbers in the z vector
z<-x+1i*y

#determine step lengths
dz<-diff(z)
#steplength
s<-Mod(dz)
data.unclass$STEPLENGTH[2:length(z)]<-s

#absolute angles
phi<-Arg(dz)
#turning angles
theta<-diff(phi)
data.unclass$TURNANGLE[3:length(z)]<-theta
assign(paste("LT",substr(tag[i],start=4,stop=5),sep=""),data.unclass)
}

fishlist2<-list(LT31,LT32,LT33,LT34,LT35,LT37,LT38,LT40,LT41,LT43,LT44)
rm(LT31,LT32,LT33,LT34,LT35,LT37,LT38,LT40,LT41,LT43,LT44)


daily<-NULL
for(i in 1:length(fishlist2)){
  fishi<-fishlist2[[i]]
  fishi<-fishi[3:dim(fishi)[1],]
  step<-aggregate(fishi$STEPLENGTH,list(fishi$DATE),sum)
  pers<-aggregate(fishi$TURNANGLE,list(fishi$DATE),function(x) mean(cos(x)))
  accel<-aggregate(fishi$ACCEL,list(fishi$DATE),function(x) mean(x,na.rm=TRUE))
  daily1<-data.frame(TRANSMITTER=rep(fishi$TRANSMITTER[1],dim(step)[1]),
                     SEX=rep(fishi$SEX[1],dim(step)[1]),
                     DATE=step[,1],
                     STEPLENGTH=step[,2],
                     PERSINDEX=pers[,2],
                     ACCEL=accel[,2]
                     )
  daily<-rbind(daily,daily1)
}

daily$DATE<-as.POSIXct(daily$DATE)
daily<-daily[daily$DATE!="2013-09-24" & daily$DATE!="2013-09-27",]
 

f<-ltspawnclust3[ltspawnclust3$SEX=="F",]
f_mean<-mean(tapply(f$SITE,f$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)
f_sd<-sd(tapply(f$SITE,f$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)


f_clust<-mean(tapply(f$CLUSTER,f$TRANSMITTER, length),na.rm=TRUE)
f_sdclust<-sd(tapply(f$CLUSTER,f$TRANSMITTER, length),na.rm=TRUE)


fa_mean<-mean(tapply(f$TOTALTIME,f$TRANSMITTER,mean),na.rm=TRUE)
fa_sd<-sd(tapply(f$TOTALTIME,f$TRANSMITTER,mean),na.rm=TRUE)

fa_sum<-mean(tapply(f$TOTALTIME,f$TRANSMITTER,sum),na.rm=TRUE)
fa_sum_sd<-sd(tapply(f$TOTALTIME,f$TRANSMITTER,sum),na.rm=TRUE)

fno<-f[f$SPAWN==0,]

fn_mean<-mean(tapply(fno$SITE,fno$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)
fn_sd<-sd(tapply(fno$SITE,fno$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)

fn_clust<-mean(tapply(fno$CLUSTER,fno$TRANSMITTER, length),na.rm=TRUE)
fn_sdclust<-sd(tapply(fno$CLUSTER,fno$TRANSMITTER, length),na.rm=TRUE)

ft_mean<-mean(tapply(fno$TOTALTIME,fno$TRANSMITTER,mean),na.rm=TRUE)
ft_sd<-sd(tapply(fno$TOTALTIME,fno$TRANSMITTER,mean),na.rm=TRUE)

ft_sum<-mean(tapply(fno$TOTALTIME,fno$TRANSMITTER,sum),na.rm=TRUE)
ft_sum_sd<-sd(tapply(fno$TOTALTIME,fno$TRANSMITTER,sum),na.rm=TRUE)

fs<-f[f$SPAWN==1,]
fts_mean<-mean(tapply(fs$TOTALTIME,fs$TRANSMITTER,mean),na.rm=TRUE)
fts_sd<-sd(tapply(fs$TOTALTIME,fs$TRANSMITTER,mean),na.rm=TRUE)

fts_sum<-mean(tapply(fs$TOTALTIME,fs$TRANSMITTER,sum),na.rm=TRUE)
fts_sum_sd<-sd(tapply(fs$TOTALTIME,fs$TRANSMITTER,sum),na.rm=TRUE)

fs_mean<-mean(tapply(fs$SITE,fs$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)
fs_sd<-sd(tapply(fs$SITE,fs$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)

fs_clust<-mean(tapply(fs$CLUSTER,fs$TRANSMITTER, length),na.rm=TRUE)
fs_sdclust<-sd(tapply(fs$CLUSTER,fs$TRANSMITTER, length),na.rm=TRUE)

m<-ltspawnclust3[ltspawnclust3$SEX=="M",]
m_mean<-mean(tapply(m$SITE,m$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)
m_sd<-sd(tapply(m$SITE,m$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)

m_clust<-mean(tapply(m$CLUSTER,m$TRANSMITTER, length),na.rm=TRUE)
m_sdclust<-sd(tapply(m$CLUSTER,m$TRANSMITTER, length),na.rm=TRUE)

ma_mean<-mean(tapply(m$TOTALTIME,m$TRANSMITTER,mean),na.rm=TRUE)
ma_sd<-sd(tapply(m$TOTALTIME,m$TRANSMITTER,mean),na.rm=TRUE)

ma_sum<-mean(tapply(m$TOTALTIME,m$TRANSMITTER,sum),na.rm=TRUE)
ma_sum_sd<-sd(tapply(m$TOTALTIME,m$TRANSMITTER,sum),na.rm=TRUE)


mno<-m[m$SPAWN==0,]
mt_mean<-mean(tapply(mno$TOTALTIME,mno$TRANSMITTER,mean),na.rm=TRUE)
mt_sd<-sd(tapply(mno$TOTALTIME,mno$TRANSMITTER,mean),na.rm=TRUE)

mt_sum<-mean(tapply(mno$TOTALTIME,mno$TRANSMITTER,sum),na.rm=TRUE)
mt_sum_sd<-sd(tapply(mno$TOTALTIME,mno$TRANSMITTER,sum),na.rm=TRUE)

mn_mean<-mean(tapply(mno$SITE,mno$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)
mn_sd<-sd(tapply(mno$SITE,mno$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)

mn_clust<-mean(tapply(mno$CLUSTER,mno$TRANSMITTER, length),na.rm=TRUE)
mn_sdclust<-sd(tapply(mno$CLUSTER,mno$TRANSMITTER, length),na.rm=TRUE)

ms<-m[m$SPAWN==1,]
mts_mean<-mean(tapply(ms$TOTALTIME,ms$TRANSMITTER,mean),na.rm=TRUE)
mts_sd<-sd(tapply(ms$TOTALTIME,ms$TRANSMITTER,mean),na.rm=TRUE)

mts_sum<-mean(tapply(ms$TOTALTIME,ms$TRANSMITTER,sum),na.rm=TRUE)
mts_sum_sd<-sd(tapply(ms$TOTALTIME,ms$TRANSMITTER,sum),na.rm=TRUE)

ms_mean<-mean(tapply(ms$SITE,ms$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)
ms_sd<-sd(tapply(ms$SITE,ms$TRANSMITTER,function(x) length(unique(x))),na.rm=TRUE)

ms_clust<-mean(tapply(ms$CLUSTER,ms$TRANSMITTER, length),na.rm=TRUE)
ms_sdclust<-sd(tapply(ms$CLUSTER,ms$TRANSMITTER, length),na.rm=TRUE)

site<-data.frame(SEX=c("Male","Female"),
                    MEANCLUST=c(m_clust,f_clust),
                    MCSD=c(m_sdclust,f_sdclust),   
                    SITEVISIT=c(m_mean,f_mean),
                    SVSD=c(m_sd,f_sd),
                    TIMEMEAN=c(ma_mean,fa_mean),
                    TIMESD=c(ma_sd,fa_sd),
                    TOTALTIME=c(ma_sum,fa_sum),
                    TOTALTIMESD=c(ma_sum_sd,fa_sum_sd),
                    
                    NOMEANCLUST=c(mn_clust,fn_clust),
                    NCSD=c(mn_sdclust,fn_sdclust),  
                    NOSITEVISIT=c(mn_mean,fn_mean),
                    NSSD=c(mn_sd,fn_sd),
                    NOTIMEMEAN=c(mt_mean,ft_mean),
                    NOTIMESD=c(mt_sd,ft_sd),
                    NOTOTALTIME=c(mt_sum,ft_sum),
                    NOTOTALTIMESD=c(mt_sum_sd,ft_sum_sd),
                    
                    MEANSPAWN=c(ms_clust,fs_clust),
                    MSSD=c(ms_sdclust,fs_sdclust),  
                    SITESPAWN=c(ms_mean,fs_mean),
                    SSSD=c(ms_sd,fs_sd),
                    SPAWNTIMEMEAN=c(mts_mean,fts_mean),
                    SPAWNTIMESD=c(mts_sd,fts_sd),
                    SPAWNTOTALTIME=c(mts_sum,fts_sum),
                    SPAWNTOTALTIMESD=c(mts_sum_sd,fts_sum_sd)
)

require(adehabitatHR)
mcpclust=NULL
ctags<-unique(ltspawnclust3$TRANSMITTER)
for(i in 1:length(fishlist)){
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))
  if(fishi$TRANSMITTER[1]%in%ctags){
    clust<-ltspawnclust3[ltspawnclust3$TRANSMITTER==fishi$TRANSMITTER[1],]
    clusters<-unique(clust$CLUSTER)
    new.xy<-xy[fishi$CLUSTER%in%clusters,]
    cp=mcp(new.xy,percent=100,unout="km2")
    #spawn clusters
    clust<-ltspawnclust3[ltspawnclust3$TRANSMITTER==fishi$TRANSMITTER[1] & ltspawnclust3$SPAWN==1,]
    clusters<-unique(clust$CLUSTER)
    new.xy<-xy[fishi$CLUSTER%in%clusters,]
    cps=mcp(new.xy,percent=100,unout="km2")
    clust<-ltspawnclust3[ltspawnclust3$TRANSMITTER==fishi$TRANSMITTER[1] & ltspawnclust3$SPAWN==0,]
    clusters<-unique(clust$CLUSTER)
    new.xy<-xy[fishi$CLUSTER%in%clusters,]
    cpn=mcp(new.xy,percent=100,unout="km2")
    mcpclust1<-data.frame(TRANSMITTER=fishi$TRANSMITTER[1],
                          SEX=fishi$SEX[1],
                          MCP=cp$area,
                          SPAWNMCP=cps$area,
                          NOMCP=cpn$area,
                          PERC=cps$area/cp$area*100)
    mcpclust<-rbind(mcpclust,mcpclust1)
  }
}
require(adehabitatHR)
xy<-SpatialPoints(cbind(spawn$CLUST.X,spawn$CLUST.Y))
all_mcp<-mcp(xy,percent=100,unout="km2")


m_area<-mean(mcpclust$MCP[mcpclust$SEX=="M"])
m_sdarea<-sd(mcpclust$MCP[mcpclust$SEX=="M"])

f_area<-mean(mcpclust$MCP[mcpclust$SEX=="F"])
f_sdarea<-sd(mcpclust$MCP[mcpclust$SEX=="F"])

mn_area<-mean(mcpclust$NOMCP[mcpclust$SEX=="M"])
mn_sdarea<-sd(mcpclust$NOMCP[mcpclust$SEX=="M"])

fn_area<-mean(mcpclust$NOMCP[mcpclust$SEX=="F"])
fn_sdarea<-sd(mcpclust$NOMCP[mcpclust$SEX=="F"])

ms_area<-mean(mcpclust$SPAWNMCP[mcpclust$SEX=="M"])
ms_sdarea<-sd(mcpclust$SPAWNMCP[mcpclust$SEX=="M"])

fs_area<-mean(mcpclust$SPAWNMCP[mcpclust$SEX=="F"])
fs_sdarea<-sd(mcpclust$SPAWNMCP[mcpclust$SEX=="F"])

m_parea<-mean(mcpclust$PERC[mcpclust$SEX=="M"])
m_sdparea<-sd(mcpclust$PERC[mcpclust$SEX=="M"])

f_parea<-mean(mcpclust$PERC[mcpclust$SEX=="F"])
f_sdparea<-sd(mcpclust$PERC[mcpclust$SEX=="F"])

#distance between spawning shoals
spawn<-NULL
ctags<-unique(ltspawnclust3$TRANSMITTER)
for(i in 1:length(ctags)){
  fishi<-ltspawnclust3[ltspawnclust3$TRANSMITTER==ctags[i],]
  t<-fishi$STARTTIME
  x<-fishi$CLUST.X
  y<-fishi$CLUST.Y
  #add x and y data as real and imaginary numbers in the z vector
  z<-x+1i*y
  #determine step lengths
  fishi$DISTCLUST[2:length(z)]<-Mod(diff(z))
  fishi$TIMECLUST[2:length(z)]<-difftime(t[c(2:length(t))],t[c(1:length(t)-1)],units="hours")
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


m_dist<-mean(spawn$DISTCLUST[spawn$SEX=="M"],na.rm=TRUE)
m_sddist<-sd(spawn$DISTCLUST[spawn$SEX=="M"],na.rm=TRUE)

f_dist<-mean(spawn$DISTCLUST[spawn$SEX=="F"],na.rm=TRUE)
f_sddist<-sd(spawn$DISTCLUST[spawn$SEX=="F"],na.rm=TRUE)

m_time<-mean(spawn$TIMECLUST[spawn$SEX=="M"],na.rm=TRUE)
m_sdtime<-sd(spawn$TIMECLUST[spawn$SEX=="M"],na.rm=TRUE)

f_time<-mean(spawn$TIMECLUST[spawn$SEX=="F"],na.rm=TRUE)
f_sdtime<-sd(spawn$TIMECLUST[spawn$SEX=="F"],na.rm=TRUE)


mn_dist<-mean(spawn$DISTNOCLUST[spawn$SEX=="M"],na.rm=TRUE)
mn_sddist<-sd(spawn$DISTNOCLUST[spawn$SEX=="M"],na.rm=TRUE)

fn_dist<-mean(spawn$DISTNOCLUST[spawn$SEX=="F"],na.rm=TRUE)
fn_sddist<-sd(spawn$DISTNOCLUST[spawn$SEX=="F"],na.rm=TRUE)

mn_time<-mean(spawn$TIMENOCLUST[spawn$SEX=="M"],na.rm=TRUE)
mn_sdtime<-sd(spawn$TIMENOCLUST[spawn$SEX=="M"],na.rm=TRUE)

fn_time<-mean(spawn$TIMENOCLUST[spawn$SEX=="F"],na.rm=TRUE)
fn_sdtime<-sd(spawn$TIMENOCLUST[spawn$SEX=="F"],na.rm=TRUE)

ms_dist<-mean(spawn$DISTSPAWNCLUST[spawn$SEX=="M"],na.rm=TRUE)
ms_sddist<-sd(spawn$DISTSPAWNCLUST[spawn$SEX=="M"],na.rm=TRUE)

fs_dist<-mean(spawn$DISTSPAWNCLUST[spawn$SEX=="F"],na.rm=TRUE)
fs_sddist<-sd(spawn$DISTSPAWNCLUST[spawn$SEX=="F"],na.rm=TRUE)

ms_time<-mean(spawn$TIMESPAWNCLUST[spawn$SEX=="M"],na.rm=TRUE)
ms_sdtime<-sd(spawn$TIMESPAWNCLUST[spawn$SEX=="M"],na.rm=TRUE)

fs_time<-mean(spawn$TIMESPAWNCLUST[spawn$SEX=="F"],na.rm=TRUE)
fs_sdtime<-sd(spawn$TIMESPAWNCLUST[spawn$SEX=="F"],na.rm=TRUE)


m_prior<-mean(spawn[match(unique(spawn$TRANSMITTER[spawn$SEX=="M"]),spawn$TRANSMITTER),]$CLUSTPRIOR)

m_sdprior<-sd(spawn[match(unique(spawn$TRANSMITTER[spawn$SEX=="M"]),spawn$TRANSMITTER),]$CLUSTPRIOR)

f_prior<-mean(spawn[match(unique(spawn$TRANSMITTER[spawn$SEX=="F"]),spawn$TRANSMITTER),]$CLUSTPRIOR)

f_sdprior<-sd(spawn[match(unique(spawn$TRANSMITTER[spawn$SEX=="F"]),spawn$TRANSMITTER),]$CLUSTPRIOR)


spawntable<-data.frame(SEX=c("Male","Female"),
                       n=c(length(unique(daily$TRANSMITTER[daily$SEX=="M"])), length(unique(daily$TRANSMITTER[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"]))),
                       DAILY_DISPLACE=c(mean(daily$STEPLENGTH[daily$SEX=="M"]), mean(daily$STEPLENGTH[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])),
                       SD_DAILY_DISPLACE=c(sd(daily$STEPLENGTH[daily$SEX=="M"]), sd(daily$STEPLENGTH[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])),
                       PI=c(mean(daily$PERSINDEX[daily$SEX=="M"]),mean(daily$PERSINDEX[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])),
                       SD_PI=c(sd(daily$PERSINDEX[daily$SEX=="M"]),sd(daily$PERSINDEX[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])),
                       ACCEL=c(mean(daily$ACCEL[daily$SEX=="M"]),mean(daily$ACCEL[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])),
                       SD_ACCEL=c(sd(daily$ACCEL[daily$SEX=="M"]),sd(daily$ACCEL[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])),
                       TOTALCLUST=c(table(ltspawnclust3$SEX)[2],table(ltspawnclust3$SEX)[1]),
                       TOTALSPAWN=c(table(ltspawnclust3$SEX[ltspawnclust3$SPAWN==1])[2],table(ltspawnclust3$SEX[ltspawnclust3$SPAWN==1])[1]),
                       
                       MEANCLUST=c(m_clust,f_clust),
                       MCSD=c(m_sdclust,f_sdclust),   
                       SITEVISIT=c(m_mean,f_mean),
                       SVSD=c(m_sd,f_sd),
                       TIMEMEAN=c(ma_mean,fa_mean),
                       TIMESD=c(ma_sd,fa_sd),
                       TOTALTIME=c(mt_sum,ft_sum),
                       TOTALTIMESD=c(mt_sum_sd,ft_sum_sd),
                       DISTBETW=c(m_dist,f_dist),
                       DSD=c(m_sddist,f_sddist),
                       TIMEBETW=c(m_time,f_time),
                       TSD=c(m_sddist,f_sddist),
                       AREA=c(m_area,f_area),
                       ASD=c(m_sdarea,f_sdarea),
                       
                       NOMEANCLUST=c(mn_clust,fn_clust),
                       NCSD=c(mn_sdclust,fn_sdclust),  
                       NOSITEVISIT=c(mn_mean,fn_mean),
                       NSSD=c(mn_sd,fn_sd),
                       NOTIMEMEAN=c(mt_mean,ft_mean),
                       NOTIMESD=c(mt_sd,ft_sd),
                       NOTOTALTIME=c(mt_sum,ft_sum),
                       NOTOTALTIMESD=c(mt_sum_sd,ft_sum_sd),
                       NODISTBETW=c(mn_dist,fn_dist),
                       NDSD=c(mn_sddist,fn_sddist),
                       NOTIMEBETW=c(mn_time,fn_time),
                       NTSD=c(mn_sddist,fn_sddist),
                       NOAREA=c(mn_area,fn_area),
                       NASD=c(mn_sdarea,fn_sdarea),
                       
                       MEANSPAWN=c(ms_clust,fs_clust),
                       MSSD=c(ms_sdclust,fs_sdclust),  
                       SITESPAWN=c(ms_mean,fs_mean),
                       SSSD=c(ms_sd,fs_sd),
                       SPAWNTIMEMEAN=c(mts_mean,fts_mean),
                       SPAWNTIMESD=c(mts_sd,fts_sd),
                       SPAWNTOTALTIME=c(mts_sum,fts_sum),
                       SPAWNTOTALTIMESD=c(mts_sum_sd,fts_sum_sd),
                       SPAWNDISTBETW=c(ms_dist,fs_dist),
                       SDSD=c(ms_sddist,fs_sddist),
                       SPAWNTIMEBETW=c(ms_time,fs_time),
                       STSD=c(ms_sdtime,fs_sdtime),
                       CLUSTPRIOR=c(m_prior,f_prior),
                       CPSD=c(m_sdprior,f_sdprior),
                       SPAWNAREA=c(ms_area,fs_area),
                       SASD=c(ms_sdarea,fs_sdarea),
                       PERCAREA=c(m_parea,f_parea),
                       PASD=c(m_sdparea,f_sdparea)
                       )


save(ltspawnclust3,temp,spawn.num,dateseq,suntimes,fishlist,daily,spawntable,spawn,mcpclust,file="SpawnClust.RData")





dynamitePlot <- function(height, error, names = NA, significance = NA, ylim = c(0,maxLim), ...){
  maxLim <- 1.1* max(mapply(sum, height, error))
  bp <- barplot(height, names.arg = names, ylim = ylim, ...)
  arrows(x0 = bp, y0 = height, y1 = height + error, angle = 90, length=0.1)
  text(x = bp, y = 3 + height + error, labels = significance, cex = 2)
}

#Spawning shoal clusters

dynamitePlot(height=site$SITEVISIT, 
             error=site$SVSD, 
             names=site$SEX,
             #significance=c("*",""),
             width = 1,
             space = 0,
             col = c("lightblue","pink"),
             xaxt="n", # suppress the default x axis
             yaxt="n", # suppress the default y axis
             ylab="Mean Spawning Sites",
             ylim=c(0,15),
             xlim=c(0,2),
             main="All Clusters")

axis(2,
     las=1)
axis(
  1, # puts the axis at the bottom
  at=c(0.5,1.5),
  labels=c("Males","Females"), # labels for each species
  lwd=0, # width of the long axis line is zero, makes invisible
  lwd.ticks=0, # width of the etick lines also zero, makes them invisible
  cex.axis=1, # offset from the axis of the labels
  mgp=c(0,0.4,0), # middle zero controls distance of labels from axis
  pos=0
)





#Total Mean Time No Spawn
dynamitePlot(height=site$NOTOTALTIME, 
             error=site$NOTOTALTIMESD, 
             names=site$SEX,
             #significance=c("*",""),
             width = 1,
             space = 0,
             col = c("lightblue","pink"),
             xaxt="n", # suppress the default x axis
             yaxt="n", # suppress the default y axis
             ylab="Total Time of Clusters (Mins)",
             ylim=c(0,7500),
             xlim=c(0,2),
             main="Total Time No Spawn Clusters")

axis(2,
     at=seq(0,7500,by=2500),
     las=1)
axis(
  1, # puts the axis at the bottom
  at=c(0.5,1.5),
  labels=c("Males","Females"), # labels for each species
  lwd=0, # width of the long axis line is zero, makes invisible
  lwd.ticks=0, # width of the etick lines also zero, makes them invisible
  cex.axis=1, # offset from the axis of the labels
  mgp=c(0,0.4,0), # middle zero controls distance of labels from axis
  pos=0
)






barplot(tapply(m$TOTALTIME,m$TRANSMITTER,sum))

#all potential spawn clusters by male
table(ltspawnclust2$TRANSMITTER[ltspawnclust2$TRANSMITTER %in% male])
#more likely spawn clusters by male
table(ltspawnclust3$TRANSMITTER[ltspawnclust3$TRANSMITTER %in% male])
table(ltspawnclust4$TRANSMITTER[ltspawnclust4$TRANSMITTER %in% male])

ltspawnclust3$TRANSMITTER

plot(shore_outline)
title("Female")
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-44",],col=4,pch=19)
f<-ltspawnclust4[ltspawnclust4$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-44",],col=2,pch=19)
legend(650500,6954500,legend=c("Low Activity","High Activity"), col=c(4,2),pch=19)



plot(shore_outline)
title("Female")
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-35",],col=4,pch=19)
f<-ltspawnclust4[ltspawnclust4$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-35",],col=2,pch=19)

plot(shore_outline)
title("Female")
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-43",],col=4,pch=19)
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-38",],col=4,pch=19)
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-41",],col=4,pch=19)


plot(shore_outline)
title("Female")
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$TRANSMITTER=="LT-38",],col=4,pch=19)







plot(shore_outline)
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f,col=rgb(1,0,0,0.5),pch=19)
m<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% male,]
points(CLUST.Y~CLUST.X,data=m,col=rgb(0,0,1,0.5),pch=19)

plot(shore_outline)
f<-ltspawnclust4[ltspawnclust4$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f,col=rgb(1,0,0,1),pch=19)
m<-ltspawnclust4[ltspawnclust4$TRANSMITTER %in% male,]
points(CLUST.Y~CLUST.X,data=m,col=rgb(0,0,1,1),pch=19,cex=0.5)



plot(shore_outline)
title("Female")
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$SPAWN==0,],col=4,pch=19)
points(CLUST.Y~CLUST.X,data=f[f$SPAWN==1,],col=2,cex=0.5,pch=19)
legend(650500,6954500,legend=c("Low Activity","High Activity"), col=c(4,2),pch=19)

plot(shore_outline)
title("Male")
m<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% male,]
points(CLUST.Y~CLUST.X,data=m[m$SPAWN==0,],col=4,pch=19)
points(CLUST.Y~CLUST.X,data=m[m$SPAWN==1,],col=2,cex=0.5,pch=19)
legend(650500,6954500,legend=c("Low Activity","High Activity"), col=c(4,2),pch=19)


plot(shore_outline)
title("Low Activity")
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$SPAWN==0,],col=rgb(1,0,0,0.4),pch=19)
m<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% male,]
points(CLUST.Y~CLUST.X,data=m[m$SPAWN==0,],col=rgb(0,0,1,0.4),pch=19)
legend(650500,6954500,legend=c("Male","Female"), col=c(4,2),pch=19)


plot(shore_outline)
title("High Activity")
m<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% male,]
points(CLUST.Y~CLUST.X,data=m[m$SPAWN==1,],col=rgb(0,0,1,0.4),pch=19)
f<-ltspawnclust3[ltspawnclust3$TRANSMITTER %in% female,]
points(CLUST.Y~CLUST.X,data=f[f$SPAWN==1,],col=rgb(1,0,0,0.4),pch=19)

legend(650500,6954500,legend=c("Male","Female"), col=c(4,2),pch=19)

plot(shore_outline)
#title("High Activity")
for (i in 1:length(ctags)){
  m<-ltspawnclust3[ltspawnclust3$TRANSMITTER==ctags[i],]
  points(CLUST.Y~CLUST.X,data=m,col=rgb(0,0,1,0.4),pch=19)
  lines(CLUST.Y~CLUST.X,data=m,col=rgb(1,0,0,0.4),lwd=4)
  
}



for (i in 1:length(ctags)){
  #plot(shore_outline)
  #title(ctags[i])
  m<-ltspawnclust3[ltspawnclust3$TRANSMITTER==ctags[i],]
  points(CLUST.Y~CLUST.X,data=m,col=rgb(0,0,1,0.4),pch=19)
  lines(CLUST.Y~CLUST.X,data=m,col=rgb(1,0,0,0.4),lwd=4)
  points(CLUST.Y~CLUST.X,data=m[m$SPAWN==1,],col=rgb(1,0,1,1))
}


plot(spawn$STARTTIME,spawn$TRANSMITTER,pch=19,cex=0.5,col=spawn$SEX,
     xlab="",
     ylab="",
     axes=FALSE
     )
box()

axis.POSIXct(1,spawn$STARTTIME)
axis(2,at=c(1:11),labels=tag,las=1)


points(CLUST.Y~CLUST.X,data=m[m$TRANSMITTER=="LT-37",],col=2)
points(CLUST.Y~CLUST.X,data=m[m$TRANSMITTER=="LT-34",],col=3)
points(CLUST.Y~CLUST.X,data=m[m$TRANSMITTER=="LT-33",],col=4)
points(CLUST.Y~CLUST.X,data=m[m$TRANSMITTER=="LT-31",],col=5)
points(CLUST.Y~CLUST.X,data=m[m$TRANSMITTER=="LT-40",],col=6)
#####clusters should be within 4 meters of bottom

bottomd<-extract(r,LT43[,c("x","y")])

#extract from bath raster contour depth - depth

#########plot clusters

x<-fishi$x
y<-fishi$y
t<-fishi$DATETIME

xyz<-cbind(x,y,t)
plot3d(x,y,t, xlim=range(x), ylim=range(y), zlim=range(t), ticktype="detailed", xlab="longitude", ylab="latitude", zlab="Time", type="l",col="grey",alpha=0.2)
for(i in 1:max(clust$cluster)){
  if(median(fishi$SHOREDIST[fishi$CLUSTER==i],na.rm=T)<40 & median(fishi$DEPTH[fishi$CLUSTER==i],na.rm=T)<5){
    new.xyz<-xyz[clust$cluster==i,]
    ch=t(convhulln(new.xyz))
    rgl.triangles(new.xyz[ch,1],new.xyz[ch,2],new.xyz[ch,3], col ="green",alpha=1)
    
  }
}

rgl.postscript("LT44_3dclust2.pdf","pdf")

rgl.open()
x.bath<-as.numeric(rownames(bath_lowres))
y.bath<-as.numeric(colnames(bath_lowres))
z.bath<-bath_lowres*10
x<-sim$UTM.X
y<-sim$UTM.Y
z<-sim$DEPTH*-10
z<-(unclass(t)-unclass(t[1]))/1000
z<-seq(1,1000,along.with=t)
rgl.surface(x.bath,y.bath,z.bath)
rgl.surface(x.bath,y.bath,matrix(0,nrow(z.bath),ncol(z.bath)),col="blue")
rgl.points(x,z,y,col="white",size=5)
plot3d(x,z,y, xlim=range(x), ylim=range(z), zlim=range(y), ticktype="detailed", xlab="longitude", ylab="time", zlab="Latitude", type="l",col="white")
rgl.lines(x,z,y)
rgl.pop()


plot3d(lon,lat,timedate, xlim=range(lon), ylim=range(lat), zlim=range(timedate), ticktype=”detailed”, xlab=”longitude”, ylab=”latitude”, zlab=”Date”, col= as.POSIXlt(daten[,"Date"])$mday, type=”l”, main=plottitle)


fishi=LT43[LT43$PART==4,]
xy<-SpatialPoints(cbind(fishi$x,fishi$y))
plot(shore_outline)
title(tag[i])
#legend(646500,6954500,legend=c("Aug 16 - 31","Sept 01 - 15", "Sept 16 - 30","Oct 01 - 15"), title="Cluster Dates", col=2:5,pch=19,cex=1.5)


lines(fishi$x,fishi$y,col=rgb(0,0,0,0.2),pch=19,cex=0.2)
points(xy,col=rgb(0,0,0,0.2),pch=19,cex=0.02)

oldspawnclust4<-ltspawnclust4
oldspawnclust3<-ltspawnclust3

###plot arc between sites

#work in progress


plot(shore_outline)
points(CLUST.Y~CLUST.X,data=spawn,pch=19,col=rgb(0,0,0,0.1))
lines(CLUST.Y~CLUST.X,data=spawn,lwd=3,col=rgb(0,0,0,0.1))

require(plotrix)

cen.x<-(spawn$CLUST.X[2]+spawn$CLUST.X[3])/2
cen.y<-(spawn$CLUST.Y[2]+spawn$CLUST.Y[3])/2
r<-sqrt((spawn$CLUST.X[2]-spawn$CLUST.X[3])^2+(spawn$CLUST.Y[2]-spawn$CLUST.Y[3])^2)/2

plot(shore_outline)
points(CLUST.Y~CLUST.X,data=spawn[2,],pch=19,col=rgb(0,0,0,0.5))
points(CLUST.Y~CLUST.X,data=spawn[3,],pch=19,col=rgb(0,0,0,0.5))
points(cen.x,cen.y,pch=19,col=rgb(1,0,0,0.5))

draw.arc(cen.x,cen.y,r,deg1=120,deg2=320)
lines(CLUST.Y~CLUST.X,data=spawn,col=rgb(0,0,0,0.5))


###draw dist and number of spawn sites on same plot



m<-spawn[spawn$SEX=="M" & spawn$SPAWN==1,]
m_sites<-as.numeric(tapply(m$SITE,m$TRANSMITTER,function(x) length(unique(x))))

f<-spawn[spawn$SEX=="F" & spawn$SPAWN==1,]
f_sites<-as.numeric(tapply(f$SITE,f$TRANSMITTER,function(x) length(unique(x))))


sites_data<-data.frame(SEX=c(rep("M",11),rep("F",11)),
                       SITES=c(m_sites,f_sites))


sites_data<-sites_data[!is.na(sites_data$SITES),]

plot(0,0,ylim=c(0,14),xlim=c(0.5,2.5),type="n",
     axes=FALSE,
     xlab="Sex",
     ylab="Number of Spawn Sites")
points(jitter(as.numeric(sites_data$SEX),factor=0.5),sites_data$SITES,pch=19,col=rgb(0,0,0,0.3))
points(spawntable$SEX,spawntable$SITESPAWN,pch=19,cex=2,col=c(rgb(0,0,1,1),rgb(1,0,0,1)))
 arrows(as.numeric(spawntable$SEX),(spawntable$SITESPAWN-spawntable$SSSD),
       as.numeric(spawntable$SEX),(spawntable$SITESPAWN+spawntable$SSSD),
       length=0,lwd=2,col=c(rgb(0,0,1,1),rgb(1,0,0,1)))
axis(1,
     at=c(1,2),
     labels=c("Females","Males"))
axis(2,at=seq(0,14,2),las=1)
box(lwd=2)


###dist between spawn events plot
dist_data<-data.frame(SEX=spawn$SEX,
                      DIST=spawn$DISTSPAWNCLUST)
dist_data<-dist_data[!is.na(dist_data$DIST),]


plot(0,0,ylim=c(-100,1500),xlim=c(0.5,2.5),type="n",
     axes=FALSE,
     xlab="Sex",
     ylab="Distance Between Spawn Sites (m)")
points(jitter(as.numeric(dist_data$SEX),factor=0.5),dist_data$DIST,pch=19,col=rgb(0,0,0,0.3))
points(spawntable$SEX,spawntable$SPAWNDISTBETW,pch=19,cex=2,col=c(rgb(0,0,1,1),rgb(1,0,0,1)))
arrows(as.numeric(spawntable$SEX),(spawntable$SPAWNDISTBETW-spawntable$SDSD),
       as.numeric(spawntable$SEX),(spawntable$SPAWNDISTBETW+spawntable$SDSD),
       length=0,lwd=2,col=c(rgb(0,0,1,1),rgb(1,0,0,1)))
axis(1,
     at=c(1,2),
     labels=c("Females","Males"))
axis(2,las=1)
box(lwd=2)
