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
data.test<-data %>%
  group_by(TRANSMITTER,SEX,SPAWN) %>%
  summarise(FREQ = length(CLUSTER),
            TOTALTIME = sum(TOTALTIME),
            VISITS = length(uniqiue(SITE)))

hist(log10(data.test$FREQ))
boxplot(FREQ~SEX*SPAWN,data=data.test)

require(lattice)

xyplot(FREQ ~ |TRANSMITTER, groups=SEX, type= c("p", "r"), data=data.test)

mod1<-gls(log(FREQ)~SEX*SPAWN,data=data.test)
#mixed effects model
mod2<-lme(log(FREQ)~SEX*SPAWN,random=~1|TRANSMITTER,data=data.test,method="REML")

anova(mod1,mod2)

#check for autocorrelation
acf(residuals(mod1))#correlation
mod3<-gls(log(FREQ)~SEX*SPAWN,correlation=corAR1(),data=data.test)
acf(residuals(mod3,type="normalized"))
anova(mod1,mod3)

mod3<-gls(log(FREQ)~SEX*SPAWN,correlation=corAR1(),data=data.test, method="ML")
Anova(mod3,type=3)



summary(mod3)


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

#### Total Duration

data.test<-data %>%
  group_by(TRANSMITTER,SEX,SPAWN) %>%
  summarise(TOTALDUR= length(CLUSTER))

hist(log10(data.test$FREQ))
boxplot(FREQ~SEX*SPAWN,data=data.test)

require(lattice)

xyplot(FREQ ~ |TRANSMITTER, groups=SEX, type= c("p", "r"), data=data.test)

mod1<-gls(log(FREQ)~SEX*SPAWN,data=data.test)
#mixed effects model
mod2<-lme(log(FREQ)~SEX*SPAWN,random=~1|TRANSMITTER,data=data.test,method="REML")

anova(mod1,mod2)

#check for autocorrelation
acf(residuals(mod1))#correlation
mod3<-gls(log(FREQ)~SEX*SPAWN,correlation=corAR1(),data=data.test)
acf(residuals(mod3,type="normalized"))
anova(mod1,mod3)

mod3<-gls(log(FREQ)~SEX*SPAWN,correlation=corAR1(),data=data.test, method="ML")
Anova(mod3,type=3)



summary(mod3)


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


