plot(shore_outline)
for(i in 1:length(fishlist)){
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))

  #legend(646500,6954500,legend=c("Aug 16 - 31","Sept 01 - 15", "Sept 16 - 30","Oct 01 - 15"), title="Cluster Dates", col=2:5,pch=19,cex=1.5)
  
  
  #lines(fishi$UTM.X,fishi$UTM.Y,col=rgb(0,0,0,0.2),pch=19,cex=0.2)
  #points(xy,col=rgb(0,0,0,0.2),pch=19,cex=0.02)
  
  for(ii in 1:max(fishi$CLUSTER)){
    new.xy<-xy[fishi$CLUSTER==ii]
    #clusters must have a median point distance < 40 m from shore and median depth < 5 m
    #if(median(fishi$SHOREDIST[fishi$CLUSTER==ii],na.rm=TRUE)<40 & median(fishi$DEPTH[fishi$CLUSTER==ii],na.rm=TRUE)<5){
    #if(length(new.xy)>=8){
    cp=mcp(new.xy,percent=100)
    #cp1=gBuffer(cp, width=50.0, quadsegs=5, capStyle="SQUARE",
    #       joinStyle="MITRE", mitreLimit=1.0)
    #plot(cp1,add=T,col=colour+1)
    plot(cp,add=TRUE,col=rgb(1,0,0,0.1),border=rgb(0,0,0,0.3))

  }
}

#plot(shore_outline,xlim=c(648000, 649000),ylim=c( 6952000, 6953000))
plot(shore_outline)
ctags<-unique(ltspawnclust3$TRANSMITTER)
for(i in 1:length(fishlist)){
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))
  if(fishi$TRANSMITTER[1]%in%ctags){
    clust<-ltspawnclust3[ltspawnclust3$TRANSMITTER==fishi$TRANSMITTER[1],]
    clusters<-unique(clust$CLUSTER)
    for(ii in 1:length(clusters)){
      new.xy<-xy[fishi$CLUSTER==clusters[ii],]
      cp=mcp(new.xy,percent=100)
      plot(cp,add=TRUE,col=rgb(1,0,0,0.4),border=rgb(0,0,0,0.3))
      
    }
  }
}


plot(shore_outline,xlim=c(648000, 649000),ylim=c( 6952000, 6953000))
plot(shore_outline)

ctags<-unique(ltspawnclust3$TRANSMITTER)
for(i in 1:length(fishlist)){
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))
  if(fishi$TRANSMITTER[1]%in%ctags){
    clust<-ltspawnclust3[ltspawnclust3$SPAWN==1 & ltspawnclust3$TRANSMITTER==fishi$TRANSMITTER[1],]
    clusters<-unique(clust$CLUSTER)
    for(ii in 1:length(clusters)){
      new.xy<-xy[fishi$CLUSTER==clusters[ii],]
      cp=mcp(new.xy,percent=100)
      plot(cp,add=TRUE,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
      
    }
  }
}

plot(shore_outline)

ctags<-unique(ltspawnclust3$TRANSMITTER[ltspawnclust3$SEX=="F"])
for(i in 1:length(fishlist)){
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))
  if(fishi$TRANSMITTER[1]%in%ctags){
    clust<-ltspawnclust3[ltspawnclust3$TRANSMITTER==fishi$TRANSMITTER[1],]
    clusters<-unique(clust$CLUSTER)
    for(ii in 1:length(clusters)){
      new.xy<-xy[fishi$CLUSTER==clusters[ii],]
      cp=mcp(new.xy,percent=100)
      plot(cp,add=TRUE,col=rgb(1,0,0,0.5),border=rgb(1,0,0,0.5))
      
    }
  }
}

ctags<-unique(ltspawnclust3$TRANSMITTER[ltspawnclust3$SEX=="M"])
for(i in 1:length(fishlist)){
  fishi=fishlist[[i]]
  xy<-SpatialPoints(cbind(fishi$UTM.X,fishi$UTM.Y))
  if(fishi$TRANSMITTER[1]%in%ctags){
    clust<-ltspawnclust3[ltspawnclust3$TRANSMITTER==fishi$TRANSMITTER[1],]
    clusters<-unique(clust$CLUSTER)
    for(ii in 1:length(clusters)){
      new.xy<-xy[fishi$CLUSTER==clusters[ii],]
      cp=mcp(new.xy,percent=100)
      plot(cp,add=TRUE,col=rgb(0,0,1,0.3),border=rgb(0,0,1,0.4))
      
    }
  }
}