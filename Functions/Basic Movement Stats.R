#format scientific notation

load("SpawnClust.RData")

#select only males or females
data<-daily[daily$MONTH==9 & (daily$SEX=="M" | daily$SEX=="F") & daily$TRANSMITTER!="LT-32",]
data$DAY<-day(as.Date(data$DATE))
data$YEAR<-as.factor(data$YEAR)
data$TRANSMITTER<-as.factor(data$TRANSMITTER)
data$SEX<-as.factor(data$SEX)



############################

#daily displacement
require(nlme)
require(multcomp)
require(lsmeans)
require(MuMIn)
require(lmtest)


#test if random effect better than linear regression
#linear model
gls1<-gls(sqrt(STEPLENGTH)~SEX+YEAR+SEX:YEAR,data=data)
#mixed effects model
mod1<-lme(sqrt(STEPLENGTH)~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER,data=data, method="REML")

anova(gls1,mod1)#random effect better model, 
#l.ratio<-anova(gls1,mod1)[[2,"L.Ratio"]] 
#0.5*(1-pchisq(l.ratio,1))
#check for autocorrelation
acf(residuals(mod1, type="normalized"))#looks like correlation

#add AR1 correlation structure to model
mod2<-lme(sqrt(STEPLENGTH)~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data, method="REML")

#the AR1 model is a better fit then the model without the correlation structure
anova(mod1,mod2)

#re-check autocorrelation
acf(residuals(mod2, type="normalized"))#looks pretty good


#re-run model 2 with maximum likelihood (ML)
mod2<-lme(sqrt(STEPLENGTH)~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data, method="ML")

summary(mod2)


#drop SEX:YEAR 
mod3<-update(mod2,.~.-SEX:YEAR)
#check significance of SEX:YEAR
anova(mod2,mod3)# looks to be insignificant and can be dropped from full model


#drop SEX
mod4.a<-update(mod3,.~.-SEX)
#check significance of SEX
anova(mod3,mod4.a)# not significant either


#drop YEAR
mod4.b<-update(mod3,.~.-YEAR)
#Check significance of YEAR
anova(mod3,mod4.b)# significant difference

#new model
mod4<-update(mod3,.~.-SEX)

#drop YEAR
mod5<-update(mod4,.~.-YEAR)
#Check significance of YEAR
anova(mod4,mod5)# significant difference

#full model
mod.dd<-lme(sqrt(STEPLENGTH)~YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data, method="REML")
summary(mod.dd)
r.squaredGLMM(mod.dd)

#check assumptions
#heterogeneity
plot(fitted(mod.dd), residuals(mod.dd, type="normalized"), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mod.dd), residuals(mod.dd, type="normalized")),col="red")
#normality
hist(residuals(mod.dd,type="normalized"))
#independence
acf(residuals(mod.dd,type="normalized"))



####Persistence Index

hist(data$PERSINDEX)


boxplot(PERSINDEX~SEX*YEAR,data=data)

#test if random effect better than linear regression
#linear model
gls1<-gls(PERSINDEX~SEX+YEAR+SEX:YEAR,data=data)
#mixed effects model
mod1<-lme(PERSINDEX~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER,data=data, method="REML")

anova(gls1,mod1)#random effect better model, 
#l.ratio<-anova(gls1,mod1)[[2,"L.Ratio"]] 
#0.5*(1-pchisq(l.ratio,1))
#check for autocorrelation
acf(residuals(mod1, type="normalized"))#looks like correlation

#add AR1 correlation structure to model
mod2<-lme(PERSINDEX~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data, method="REML")

#the AR1 model is a better fit then the model without the correlation structure
anova(mod1,mod2)

#re-check autocorrelation
acf(residuals(mod2, type="normalized"))#looks pretty good

mod2<-lme(PERSINDEX~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data, method="ML")

#test significance of intereaction
mod2.a<-update(mod2,.~.-SEX:YEAR)# drop interaction
anova(mod2,mod2.a)#not significant remove

mod3<-mod2.a

#drop sex
mod3.a<-update(mod3,.~.-SEX)
anova(mod3,mod3.a)#significant
#drop year
mod3.b<-update(mod3,.~.-YEAR)
anova(mod3,mod3.b)#significant

#best model

mod.pi<-lme(PERSINDEX~SEX+YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data, method="REML")
summary(mod.pi)
r.squaredGLMM(mod.pi)

#check assumptions
#heterogeneity
plot(fitted(mod.pi), residuals(mod.pi, type="normalized"), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mod.pi), residuals(mod.pi, type="normalized")),col="red")
#normality
hist(residuals(mod.pi,type="normalized"))
#independence
acf(residuals(mod.pi,type="normalized"))



###Acceleration

#accelerometers not available in 2012
data.ac<-data[data$YEAR!=2012,]
data.ac$YEAR<-droplevels(data.ac$YEAR)
#remove fish without accelerometers
data.ac<-data.ac[!is.na(data.ac$ACCELMEAN),]

hist(data.ac$ACCELMEAN)
#log transform for normality
hist(log(data.ac$ACCELMEAN))
boxplot(log(ACCELMEAN)~SEX*YEAR,data.ac)
#test if random effect better than linear regression
#linear model
gls1<-gls(log(ACCELMEAN)~SEX+YEAR+SEX:YEAR,data=data.ac)
#mixed effects model
mod1<-lme(log(ACCELMEAN)~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER,data=data.ac, method="REML")

anova(gls1,mod1)#random effect better model, 
#l.ratio<-anova(gls1,mod1)[[2,"L.Ratio"]] 
#0.5*(1-pchisq(l.ratio,1))
#check for autocorrelation
acf(residuals(mod1, type="normalized"))#looks like correlation

#add AR1 correlation structure to model
mod2<-lme(log(ACCELMEAN)~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data.ac, method="REML")

#the AR1 model is a better fit then the model without the correlation structure
anova(mod1,mod2)

#re-check autocorrelation
acf(residuals(mod2, type="normalized"))#looks pretty good

mod2<-lme(log(ACCELMEAN)~SEX+YEAR+SEX:YEAR, random=~1|TRANSMITTER, correlation=corAR1() ,data=data.ac, method="ML")

#test significance of intereaction
mod2.a<-update(mod2,.~.-SEX:YEAR)# drop interaction
anova(mod2,mod2.a)#not significant remove

mod3<-mod2.a

#drop sex
mod3.a<-update(mod3,.~.-SEX)
anova(mod3,mod3.a)#not significant but on the margin
l.ratio<-anova(mod3,mod3.a)[[2,"L.Ratio"]] 
0.5*(1-pchisq(l.ratio,1))#still not significant


#drop year
mod3.b<-update(mod3,.~.-YEAR)
anova(mod3,mod3.b)#not significant and greater

mod4<-mod3.b
#drop sex
mod4.a<-update(mod4,.~.-SEX)
anova(mod4,mod4.a)#not significant

#best model

mod.ac<-lme(log(ACCELMEAN)~1, random=~1|TRANSMITTER, correlation=corAR1() ,data=data.ac, method="REML")
summary(mod.ac)
r.squaredGLMM(mod.ac)

#check assumptions
#heterogeneity
plot(fitted(mod.ac), residuals(mod.ac, type="normalized"), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mod.ac), residuals(mod.ac, type="normalized")),col="red")
#normality
hist(residuals(mod.ac,type="normalized"))
#independence
acf(residuals(mod.ac,type="normalized"))

library(car)
Anova(update(mod2,.~.-YEAR-SEX:YEAR))


