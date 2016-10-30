tags<-unique(fish$TRANSMITTER)
par(mfrow=c(2,3))
i=1
plot(UTM.Y~UTM.X,data=fish[fish$TRANSMITTER==tags[i],],main=tags[i])
plot(rep(1,length(DATETIME))~DATETIME,data=fish[fish$TRANSMITTER==tags[i],],main=tags[i])
plot(DEPTH~DATETIME,data=fish[fish$TRANSMITTER==tags[i],],main=tags[i])
#2012
plot(shore_outline, main="2012")
points(UTM.Y~UTM.X,data=fish[fish$TRANSMITTER==tags[i] & 
                               fish$DATETIME<"2013-01-01",])
#2013
plot(shore_outline, main="2013")
points(UTM.Y~UTM.X,data=fish[fish$TRANSMITTER==tags[i] & 
                               fish$DATETIME>"2013-01-01" &
                               fish$DATETIME<"2014-01-01",])
#2014
plot(shore_outline, main="2014")
points(UTM.Y~UTM.X,data=fish[fish$TRANSMITTER==tags[i] & 
                               fish$DATETIME>"2014-01-01",])
i=i+1
