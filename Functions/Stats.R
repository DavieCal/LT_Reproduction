#format scientific notation

library(lme4)
load("SpawnClust.RData")
#select only september 
data<-daily[daily$MONTH==9,]
#select only males or females
data<-daily[daily$MONTH==9 & (daily$SEX=="M" | daily$SEX=="F"),]

library(fitdistrplus)
library(car)
library(MASS)

y<-data$STEPLENGTH
descdist(y)
fit.norm<-fitdist(y,"norm")
fit.gamma<-fitdist(y,"gamma")




#Daily displacement
qqp(data$STEPLENGTH,"gaussian")
hist(sqrt(data$STEPLENGTH))


boxplot(sqrt(STEPLENGTH)~SEX*YEAR,data)

data$YEAR<-as.factor(data$YEAR)
data$TRANSMITTER<-as.factor(data$TRANSMITTER)
data$SEX<-as.factor(data$SEX)
#REML= FALSE necessary to do likelihood ratio test
fit.dd.full<-lmer(sqrt(STEPLENGTH)~SEX+(1|TRANSMITTER)+(1|YEAR),data=data, REML=FALSE)
fit.dd.null<-lmer(sqrt(STEPLENGTH)~1+(1|TRANSMITTER)+(1|YEAR),data=data, REML=FALSE)

anova(fit.dd.full,fit.dd.null)
qqnorm(resid(fit.dd.full))
summary(fit.dd.full)
#look at residuals
plot(fit.dd.full)
#normality of residuals
hist(residuals(fit.dd.full))

coef(fit.dd.full)



#Persistence index

hist(data$PERSINDEX)

boxplot(PERSINDEX~SEX*YEAR,data)

fit.pi.full<-lmer(PERSINDEX~SEX+(1|TRANSMITTER)+(1|YEAR),data=data, REML=FALSE)
fit.pi.null<-lmer(PERSINDEX~1+(1|TRANSMITTER)+(1|YEAR),data=data, REML=FALSE)
anova(fit.pi.full,fit.pi.null)
summary(fit.pi.full)
hist(residuals(fit.pi.full))

#Acceleration

data.ac<-data[data$YEAR!=2012,]
levels(data.ac$YEAR)<-c("2013","2014")
hist(log(data.ac$ACCEL))
boxplot(ACCEL~SEX*YEAR,data.ac)


fit.ac.full<-lmer(log(ACCEL)~SEX+YEAR+(1|TRANSMITTER),data=data.ac, REML=FALSE)
fit.ac.sex<-lmer(log(ACCEL)~SEX+(1|TRANSMITTER),data=data.ac, REML=FALSE)
fit.ac.year<-lmer(log(ACCEL)~YEAR+(1|TRANSMITTER),data=data.ac, REML=FALSE)
fit.ac.null<-lmer(log(ACCEL)~(1|TRANSMITTER),data=data.ac, REML=FALSE)

anova(fit.ac.full,fit.ac.sex)
anova(fit.ac.full,fit.ac.year)
anova(fit.ac.full,fit.ac.sex)
summary(fit.ac.sex)
plot(fit.ac.full)
hist(residuals(fit.ac.full))

coef(fit.ac.full)

############################




