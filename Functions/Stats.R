#format scientific notation

require(car)
spawn$SPAWN<-as.factor(spawn$SPAWN)

#total clust
clust<-as.data.frame(table(spawn$SEX,spawn$TRANSMITTER))
clust<-clust[clust$Freq!=0,-2]
m_min<-min(clust$Freq[clust$Var1=="M"])
m_max<-max(clust$Freq[clust$Var1=="M"])

f_min<-min(clust$Freq[clust$Var1=="F"])
f_max<-max(clust$Freq[clust$Var1=="F"])


#clust sex and activity
clust<-as.data.frame(table(spawn$SEX,spawn$SPAWN,spawn$TRANSMITTER))
clust<-clust[clust$Freq!=0,]

#hist(log10(clust$Freq))

#shapiro.test(log1(clust$Freq))
#leveneTest(log(clust$Freq)~clust$Var1*clust$Var2)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

clustmod<-lm(log(Freq)~Var1*Var2,data=clust)



aclust.test<-Anova(clustmod, type=3)
#total duration

totaldur<-aggregate(TOTALTIME~SEX+SPAWN+TRANSMITTER,data=spawn,sum)


#hist(log(totaldur$TOTALTIME))
#shapiro.test(log(totaldur$TOTALTIME))
require(car)
#leveneTest(log(totaldur$TOTALTIME)~totaldur$SEX*totaldur$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

tdurmod<-lm(log(TOTALTIME)~SEX*SPAWN,data=totaldur)

tdurp.test<-Anova(tdurmod, type=3)



#duration per cluster

dur<-data.frame(SEX=spawn$SEX,
                SPAWN=as.factor(spawn$SPAWN),
                TRANSMITTER=spawn$TRANSMITTER,
                DURATION=spawn$TOTALTIME)


#hist(log(dur$DURATION))
#shapiro.test(log(dur$DURATION))
#bartlett.test(log(dur$DURATION)~dur$SEX)

require(car)
#leveneTest(log(dur$DURATION)~dur$SEX*dur$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

durmod<-lm(log(DURATION)~SEX*SPAWN,data=dur)

adurp.test<-Anova(durmod, type=3)



#unique site visits
require(car)
visit<-aggregate(SITE~SEX+SPAWN+TRANSMITTER,data=spawn,function(x) length(unique(x)))

#hist(log(visit$SITE))
#shapiro.test(log(visit$SITE))
#leveneTest(log(visit$SITE)~visit$SEX*visit$SPAWN)

#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

vismod<-lm(log(SITE)~SEX*SPAWN,data=visit)

avisit.test<-Anova(vismod, type=3)


#dist between

dist<-data.frame(SEX=spawn$SEX,
                 SPAWN=as.factor(spawn$SPAWN),
                 TRANSMITTER=spawn$TRANSMITTER,
                 DISTCLUST =spawn$DISTCLUST,
                 TIMECLUST=spawn$TIMECLUST)



#hist(log(dist$DISTCLUST))

require(car)
#leveneTest((log(dist$DISTCLUST))~dist$SEX*dist$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

distmod<-lm(log(DISTCLUST)~SEX*SPAWN,data=dist)

adist.test<-Anova(distmod,type=3)


#time between

#hist(1/(dist$TIMECLUST)^(1/3))

#shapiro.test(1/(dist$TIMECLUST)^(1/3))
#qqnorm(1/((dist$TIMECLUST)^(1/3)))
require(car)
#leveneTest(log(dur$DURATION)~dur$SEX*dur$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

timemod<-lm(1/(TIMECLUST)^(1/3)~SEX*SPAWN,data=dist)

atime.test<-Anova(timemod, type=3)


#sites prior to spawn
prior<-data.frame(SEX=spawn$SEX,
                  SPAWN=as.factor(spawn$SPAWN),
                  TRANSMITTER=spawn$TRANSMITTER,
                  CLUSTPRIOR =spawn$CLUSTPRIOR)
prior<-aggregate(CLUSTPRIOR~SEX+TRANSMITTER,data=spawn,mean)


prior.test<-wilcox.test(CLUSTPRIOR~SEX,data=prior)





#MCP

smcp<-mcpclust[,-c(3,5,6)]
smcp$SPAWN<-1
names(smcp)<-c("TRANSMITTER","SEX","MCP","SPAWN")
nmcp<-mcpclust[,-c(3,4,6)]
nmcp$SPAWN<-0
names(nmcp)<-c("TRANSMITTER","SEX","MCP","SPAWN")
mcp<-rbind(smcp,nmcp)
mcp$SPAWN<-as.factor(mcp$SPAWN)

#shapiro.test(log(mcp$MCP))
require(car)
#leveneTest(log(MCP)~SEX*SPAWN,data=mcp)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))
areamod<-lm(log(MCP)~SEX*SPAWN,data=mcp)
aarea.test<-Anova(areamod,type=3)

