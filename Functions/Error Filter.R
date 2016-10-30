
#HPE Filter determine by Meckley et al 2014 method

#HPE filter to remove points with a 2DRMS > 27 m
fish<-pos[pos$HPE<106.1,]

#data retention 94%
hpe.dat.ret<-round(dim(fish)[1]/dim(pos)[1],digits=2)*100

#check out where the bad points are
plot(shore_outline)
points(UTM.Y~UTM.X,data=pos[pos$HPE>106,],col=rgb(1,0,0,0.01),pch=19)
points(UTM.Y~UTM.X,data=pos[pos$HPE<106,],col=rgb(0,0,0,0.01),pch=19)

#change accel depths na to -1
fish$DEPTH[is.na(fish$DEPTH)]=-1

####BioFilter####

source("Functions/Bio Filter Function.R")

#split dataframe into a list of tag data frames
data_unclean<-split(fish,fish$TRANSMITTER)
#data_unclean= VPS data for single tag
#min_interval= minimum ping interval
#run biolfilter on list
data.clean<-lapply(data_unclean,function(x){bioFilter(data=x, min_interval=80,velocity_cut=1.2)})
#bind the lists back into a dataframe
fish_correct<-do.call(rbind,data.clean)
fish_correct<-fish_correct[order(fish_correct$DATETIME),]
fish_correct$DEPTH[fish_correct$DEPTH==-1]=NA
fish_final<-fish_correct[fish_correct$DATASCORE==3,]

#87% data retention
bio.data.ret<-dim(fish_final)[1]/dim(pos)[1]


fish<-fish_final


save(fish,file="~/Data/filter.fish.RData")

#remove clutter
rm(data_unclean)
rm(data.clean)
rm(fish_correct)
rm(pos)
rm(fish_final)
