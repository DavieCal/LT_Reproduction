library(rgdal)
library(PBSmapping)
library(maptools)

pts<-read.csv("~/Data/GIS/spawnsites.csv")
pts<-load_coords(data=pts, date=FALSE)

for(i in 1:dim(pts)[1]){
  d<-pointDistance(pts[i,][,c("UTM.X","UTM.Y")],pts[,c("UTM.X","UTM.Y")],allpairs=TRUE,longlat=FALSE)
  if(any(d<70 & d!=0)){
    print(pts$SITE[i])
    #near sites
    print(pts$SITE[d<70 & d!=0])
  }
}

#find duplicates visually
#i=2

#plot(UTM.Y~UTM.X,data=pts[i,],col="white",
#     xlim=c(pts$UTM.X[i]-200,pts$UTM.X[i]+200),
#     ylim=c(pts$UTM.Y[i]-200,pts$UTM.Y[i]+200))
#plot(shore_outline,add=TRUE)
#points(UTM.Y~UTM.X,data=pts[1:i-1,],pch=19,cex=0.5)
#text(pts$UTM.X[1:i-1],pts$UTM.Y[1:i-1],pts$SITE[1:i-1])
#points(UTM.Y~UTM.X,data=pts[i,],col="red")
#title(pts$SITE[i-1])
#pts[i-1,]

#i=i+1

#remlist<-NULL
#remlist<-rbind(remlist,i-1)

remlist<-c(28, 29, 31, 33, 35, 38, 40, 42, 44, 48, 
51, 55, 56, 57, 64, 66, 75, 77, 81, 82)

spawn_pts<-pts[-remlist,]

#order spawn sites
spawn_pts<-spawn_pts[order(as.numeric(substr(spawn_pts$SITE,2,3))),]

#remove NA's
spawn_pts$EGGS[is.na(spawn_pts$EGGS)]<-0
coordinates(spawn_pts)<-~UTM.X+UTM.Y
proj4string(spawn_pts)<-CRS("+proj=utm +zone=11 +ellps=WGS84")

egg_pts<-spawn_pts[spawn_pts$EGGS==1,]

write.csv(egg_pts,"~/Data/GIS/eggsites.csv",row.names=FALSE)
write.csv(spawn_pts,"~/Data/GIS/spawnsites.csv",row.names=FALSE)
plot(shore_outline)
points(UTM.Y~UTM.X,data=spawn_pts,pch=19)
points(UTM.Y~UTM.X,data=spawn_pts[spawn_pts$EGGS==1,],pch=19,col="red")
