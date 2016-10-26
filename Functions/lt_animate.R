###############################
#ANIMATE SCRIPT
###############################

#Author: Dave Callaghan
#Date Created:22-04-2014
#Date Edited:22-04-2014
#
###############################
#Libraries
require(maptools)



###############################
#Set Variables

#Save current parsettings
opar <- par() 




spawn<-fishbind[fishbind$SEASON=="SPAWN",]


male<-c("LT-31","LT-33","LT-34","LT-36","LT-37","LT-40","LT-42")
female<-c("LT-32","LT-35","LT-38","LT-41","LT-43","LT-44")

spawn_male<-spawn[grep(paste(male, collapse="|"), spawn$TRANSMITTER),]
spawn_female<-spawn[grep(paste(female, collapse="|"), spawn$TRANSMITTER),]


fish1<-spawn[spawn$TRANSMITTER=="LT-44",]
#assign x and y coordinates
x<-fish1$UTM.X
y<-fish1$UTM.Y
#assign depth variable
d1<-fish1$DEPTH
#assign acceleration variable
a1<-fish1$ACCEL
#assign timestamp
t1<-fish1$DATETIME


#assign shore outline to lake variable
poly = shore_outline

#assign bathymetric raster
#r<-r


###############################
# Start processing
#assign x and y to real and imaginary values
z<-x+1i*y

z1<-x+1i*y


#b <- extract(r, cbind(x,y))
#s<-extract(shore.dist,cbind(x,y))


plot(poly)
points(fish1$UTM.X,fish1$UTM.Y,col=rgb(0,0,0,0.05),
       pch=0,cex=0.2)
#points(fish_pos$UTM.X,fish_pos$UTM.Y,col=rgb(0,0,0,0.05),pch=0,cex=0.2)



points(z,col=rgb(1,0,0,0.05),pch=0,cex=0.2)
points(fish_male$UTM.X,fish_male$UTM.Y,col=rgb(0,0,1,0.05),pch=0,cex=0.2)
points(fish_female$UTM.X,fish_female$UTM.Y,col=rgb(1,0,0,0.05),pch=0,cex=0.2)



quartz()

  
m <- rbind(c(1,1,1,1),c(1,1,1,1),c(1,1,1,1), c(2,2,3,3))
layout(m)
for(i in seq(1,length(z),by=20)){
  par(mar = c(0, 0, 3, 0))
  plot(poly,axes=F,xlab="",ylab="",asp=1)
  lines(z[1:i],type="o",col=rgb(0,0,1,.01))
  points(z[i],cex=2,pch=19,col="blue")
  lines(z1[1:i],type="o",col=rgb(1,0,0,.01))
  points(z1[i],cex=2,pch=19,col="red")
  title(t[i])
  par(mar = c(3, 3, 1, 2))
  plot(-d[1:i]~t[1:i],
       xlim=c((t[i]-(60*60*24)),(t[i]+(60*60*24))),
       ylim=c(-32,0),
       type="o",
       col="blue")
  lines(-d1[1:i]~t[1:i],col="red",type="o")
  title("Depth")
  plot(a[1:i]~t[1:i],
       xlim=c((t[i]-(60*60*24)),(t[i]+(60*60*24))),
       ylim=c(0,5),
       type="o",col="blue")
  lines(a1[1:i]~t[1:i],col="red",type="o")
  title("Acceleration")
  #Sys.sleep(0.07)
}
dev.off()

plot(-d[1:7000]~t[1:7000],type="l",ylim=c(-32,0))
plot(-d[7000:14000]~t[7000:14000],type="l",ylim=c(-32,0))

plot(-d[7000:14000]~t[7000:14000],type="l")
lines((a-20)~t,type="l",col="red")
lines(b[7000:14000]~t[7000:14000],type="l",col="green")

library(animation)
movement = function (ani.pause=F) {
  
  m <- rbind(c(1, 1,1,1),c(1,1,1,1),c(1,1,1,1), c(2,2,3,3))
  layout(m)
  for(i in seq(1,length(z),by=20)){
    par(mar = c(0, 0, 3, 0))
    plot(poly,axes=F,xlab="",ylab="",asp=1)
    lines(z[1:i],type="o",col=rgb(0,0,1,.05))
    points(z[i],cex=2,pch=19,col="blue")
    lines(z1[1:i],type="o",col=rgb(1,0,0,.05))
    points(z1[i],cex=2,pch=19,col="red")
    title(t[i])
    par(mar = c(3, 3, 1, 2))
    plot(-d[1:i]~t[1:i],
         xlim=c((t[i]-(60*60*24)),(t[i]+(60*60*24))),
         ylim=c(-32,0),
         type="o",
         col="blue")
    lines(-d1[1:i]~t[1:i],col="red",type="o")
    title("Depth")
    plot(a[1:i]~t[1:i],
         xlim=c((t[i]-(60*60*24)),(t[i]+(60*60*24))),
         ylim=c(0,5),
         type="o",col="blue")
    lines(a1[1:i]~t[1:i],col="red",type="o")
    title("Acceleration")
    #Sys.sleep(0.07)
  }
}
ani.options(interval = .15)
quartz()
movement()
dev.off()  
ani.options(ani.width=800, ani.height=800, interval=.25)
setwd("~/Documents/R-Directory/AniMove")


# Start processing
print("Set variables and start processing")
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")
saveHTML(movement( ani.pause=T), img.name = "lt", htmlfile="ltspawn.html",outdir = getwd(), replace=T)

# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")



plot(-d~t,
     #xlim=c((t[i]-5),(t[i]+5)),
     #ylim=c(0,5),
     type="o")
lines(b[1:14000]~t[1:14000],col="red")
t[1]




dz<-diff(z)
#steplength
s<-Mod(dz)

#Distribution of step lengths
hist(s, col="grey", bor="darkgrey", freq=FALSE)
lines(density(s), col=2, lwd=2)

#distribution of time steps
dt<-diff(t)
tt<-unclass(dt)
hist(tt, col="grey", bor="darkgrey", freq=FALSE)
lines(density(tt), col=2, lwd=2)

#absolute angle
phi<-Arg(dz)
#turning angle
theta<-diff(phi)

summary(s)
hist(theta)




hist(90-(Arg(dz)*180)/pi)

require(circular)
theta<-as.circular(theta)
phi<-as.circular (phi)
rose.diag(phi, bins=16)



par(mar = c(0, 0, 3, 0))
plot(poly,axes=F,xlab="",ylab="",asp=1)
lines(z[i:j],type="o",col=rgb(0,0,0,.1))
#lines(z[(i-40):i],type="o",col=rgb(0,1,0,.10))
points(z[j],cex=2,col="red")
title(t[i])
par(mar = c(3, 3, 1, 2))
plot(-d[i:j]~t[i:j],
     xlim=c(t[i],t[j]),
     #ylim=c(-31,0),
     type="o")
lines(b[i:j]~t[i:j],col="red")
title("Depth")
plot(a[i:j]~t[i:j],
     xlim=c(t[i],t[j]),
     ylim=c(0,5),
     type="o")
title("Acceleration")


