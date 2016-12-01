#format scientific notation

load("SpawnClust.RData")

clustspawn<-spawn[spawn$TRANSMITTER %in% paste0("LT-",31:44),]

#select only males or females
data<-clustspawn[clustspawn$TRANSMITTER!="LT-32",]
data$DAY<-day(as.Date(data$DATE))

data$SPAWN<-as.factor(data$SPAWN)
data$TRANSMITTER<-as.factor(data$TRANSMITTER)
data$TRANSMITTER<-droplevels(data$TRANSMITTER)
data$SEX<-as.factor(data$SEX)
data$SEX<-droplevels(data$SEX)



############################

#clusters
require(nlme)
require(MASS)
require(lme4)
require(multcomp)
require(lsmeans)
require(MuMIn)
library(dplyr)

data.test<-data %>%
  group_by(TRANSMITTER,SEX,SPAWN) %>%
  summarise(FREQ= length(CLUSTER))

hist(log10(data.test$FREQ))
boxplot(FREQ~SEX*SPAWN,data=data.test)

require(lattice)

xyplot(FREQ ~ SPAWN|TRANSMITTER, groups=SEX, type= c("p", "r"), data=data.test)


#test if random effect better than linear regression
#linear model
glm1<-glm(FREQ~SEX*SPAWN,family=poisson,data=data.test)
#mixed effects model
#set optimizer for lme
mod1<-glmer(FREQ~SEX*SPAWN + (1|TRANSMITTER), family=poisson,data=data.test)

anova(mod1,glm1)

#check for autocorrelation
acf(residuals(mod1))#looks like correlation


#re-run model 2 with maximum likelihood (ML)
mod1<-glmer(FREQ~SEX+SPAWN+ SEX:SPAWN  + (1|TRANSMITTER), family=poisson,data=data.test,REML=FALSE)
summary(mod1)


#drop SEX
mod1.a<-update(mod1,.~.-SEX)
#check significance of SEX
anova(mod1,mod1.a)# looks to be insignificant and can be dropped from full model


#drop SPAWN
mod1.b<-update(mod1,.~.-SPAWN)
#check significance of SPAWN
anova(mod1,mod1.b)# not significant either

#drop interaction
mod1.c<-update(mod1,.~.-SEX:SPAWN)
#check significance of SEX:SPAWN
anova(mod1,mod1.c)# not significant either



mod2<-mod1.b
#drop SEX
mod2.a<-update(mod2,.~.-SEX)
#check significance of SEX
anova(mod2,mod2.a)# looks to be insignificant and can be dropped from full model


#drop SEX
mod2.b<-update(mod2,.~.-SEX:SPAWN)
#check significance of SEX
anova(mod2,mod2.b)# looks to be insignificant and can be dropped from full model

mod3<-mod2.a

mod3.a<-update(mod3,.~.-SEX:SPAWN)
anova(mod3,mod3.a)

#new model
mod5<-update(mod4,.~.-SEX)

#drop YEAR
mod6<-update(mod5,.~.-YEAR)
#Check significance of YEAR
anova(mod5,mod6)# significant difference

#full model
mod.dd<-lme(sqrt(STEPLENGTH)~YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data, method="REML")
summary(mod.dd)
r.squaredGLMM(mod.dd)

#check assumptions
#heterogeneity
plot(fitted(mod1), residuals(mod1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mod1), residuals(mod1)),col="red")
#normality
hist(residuals(mod1))
#independence
acf(residuals(mod1))

####
