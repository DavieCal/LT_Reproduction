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
require(MuMIn)
require(dplyr)
require(car)
require(MASS)
require(mgcv)
data.test<-data %>%
  group_by(TRANSMITTER,SEX,SPAWN) %>%
  summarise(FREQ = length(CLUSTER),
            TOTALTIME = sum(TOTALTIME),
            VISITS = length(unique(SITE)))

hist(log(data.test$FREQ))
boxplot(FREQ~SEX*SPAWN,data=data.test)

require(lattice)

xyplot(FREQ ~ SPAWN|TRANSMITTER, groups=SEX, type= c("p", "r"), data=data.test)

mod1<-glm.nb(FREQ~SEX*SPAWN,data=data.test)

#mixed effects model
mod2<-glmer(FREQ~SEX*SPAWN+(1|TRANSMITTER),family=negative.binomial(1.3),data=data.test)
#mixed effects model

mod3<-gamm(FREQ~SEX*SPAWN, random=list(TRANSMITTER=~1),
              family=negative.binomial(1.3),
              data=data.test)

#check for autocorrelation
acf(residuals(mod3$lme))#correlation


Anova(mod2,type=3)


Anova(mod2,type=3)

plot(allEffects(mod2))


#check assumptions
#heterogeneity
plot(mod3)
plot(fitted(mod$lme), resid(mod3$lme,type="normalized"), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mod3$lme,type="normalized"), resid(mod3$lme)),col="red")
plot(data.test$SEX,residuals(mod3))
plot(data.test$SPAWN,residuals(mod3))
#normality
hist(residuals(mod3))
#independence
acf(residuals(mod3))

#### Total Duration


hist(data.test$TOTALTIME)
qqnorm(resid(mod2))
boxplot(log(TOTALTIME)~SEX*SPAWN,data=data.test)

require(lattice)

xyplot(log(TOTALTIME) ~ SPAWN|TRANSMITTER, groups=SEX, type= c("p", "r"), data=data.test)

mod1<-gls(log(TOTALTIME)~SEX*SPAWN,data=data.test)
#mixed effects model
mod2<-lme(log(TOTALTIME)~SEX*SPAWN,random=~1|TRANSMITTER,data=data.test,method="REML")

anova(mod1,mod2)

#check for autocorrelation
acf(residuals(mod2))#correlation

mod2<-lme(TOTALTIME~SEX*SPAWN,random=~1|TRANSMITTER,data=data.test,method="REML")

Anova(mod2,type=3)

mod3<-lme(TOTALTIME~SEX*SPAWN,random=~1|TRANSMITTER,
          weights=varIdent(form=~1|SEX),data=data.test,method="REML")

anova(mod2,mod3)

mod4<-lme(TOTALTIME~SEX+SPAWN,random=~1|TRANSMITTER,data=data.test,method="ML")

anova(mod2,mod3)

mod4<-lme(log(TOTALTIME)~SEX,random=~1|TRANSMITTER,data=data.test,method="ML")

anova(mod3,mod4)

mod5<-lme(log(TOTALTIME)~SPAWN,random=~1|TRANSMITTER,data=data.test,method="ML")

anova(mod3,mod5)

Anova(mod3,type=2)


#check assumptions
#heterogeneity
plot(fitted(mod3), residuals(mod3,type="normalized"), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mod3), residuals(mod3,type="normalized")),col="red")
plot(data.test$SEX,residuals(mod3,type="normalized"))
plot(data.test$SPAWN,residuals(mod3,type="normalized"))
#normality
hist(residuals(mod3,type="normalized"))
#independence
acf(residuals(mod3,type="normalized"))

plot(allEffects(mod3))

