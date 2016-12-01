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
require(MuMIn)
require(lattice)

xyplot(sqrt(STEPLENGTH) ~ DAY|TRANSMITTER, groups=YEAR, type= c("p", "r"), data=data)

#test if random effect better than linear regression
#linear model
gls1<-gls(sqrt(STEPLENGTH)~SEX+DAY+I(DAY^2)+YEAR,data=data)
#mixed effects model
mod1<-lme(sqrt(STEPLENGTH)~SEX+DAY+I(DAY^2)+YEAR, random=~1|TRANSMITTER,data=data, method="REML")

#check variance inflation factor
vif(mod1)
#centre days
data$cDAY<-data$DAY-mean(data$DAY)

#linear model
gls1<-gls(sqrt(STEPLENGTH)~SEX+cDAY+I(cDAY^2)+YEAR,data=data)
#mixed effects model
mod1<-lme(sqrt(STEPLENGTH)~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER,data=data, method="REML")
#re-check variance inflation factor
vif(mod1)

#test is random effects improves model
anova(gls1,mod1)#random effect better model,

acf(residuals(mod1, type="normalized"))#looks like correlation

#add AR1 correlation structure to model
mod2<-lme(sqrt(STEPLENGTH)~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER, correlation=corAR1(form=~1|TRANSMITTER/YEAR) ,data=data, method="REML")

#the AR1 model is a better fit then the model without the correlation structure
anova(mod1,mod2)

#re-check autocorrelation
acf(residuals(mod2, type="normalized"))#looks pretty good


#re-run model 2 with maximum likelihood (ML)
mod2<-lme(sqrt(STEPLENGTH)~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER, correlation=corAR1(form=~1|TRANSMITTER/YEAR) ,data=data, method="ML")

summary(mod2)


#drop SEX
mod2.a<-update(mod2,.~.-SEX)
#check significance of SEX
anova(mod2,mod2.a)# looks to be insignificant and can be dropped from full model


#drop DAY
mod2.b<-update(mod2,.~.-cDAY)
#check significance of DAY
anova(mod2,mod2.b)# significant!

#drop DAY^2
mod2.c<-update(mod2,.~.-I(cDAY^2))
#check significance of DAY^2
anova(mod2,mod2.c)# significant!


#drop YEAR
mod2.d<-update(mod2,.~.-YEAR)
#Check significance of YEAR
anova(mod2,mod2.d)# significant difference



#full model
mod.dd<-mod2.a
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
xyplot(PERSINDEX ~ DAY|TRANSMITTER, groups=YEAR, type= c("p", "r"), data=data)


#test if random effect better than linear regression
#linear model
gls1<-gls(PERSINDEX~SEX+cDAY+I(cDAY^2)+YEAR,data=data)
#mixed effects model
mod1<-lme(PERSINDEX~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER,data=data, method="REML")

#check variance inflation factor
vif(mod1)
anova(gls1,mod1)#random effect better model,
#l.ratio<-anova(gls1,mod1)[[2,"L.Ratio"]] 
#0.5*(1-pchisq(l.ratio,1))
#check for autocorrelation
acf(residuals(mod1, type="normalized"))#looks like correlation

#add AR1 correlation structure to model
mod2<-lme(PERSINDEX~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER, correlation=corAR1(form=~1|TRANSMITTER/YEAR) ,data=data, method="REML")


#the AR1 model is a better fit then the model without the correlation structure
anova(mod1,mod2)

#re-check autocorrelation
acf(residuals(mod2, type="normalized"))#looks pretty good

mod2<-lme(PERSINDEX~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER, correlation=corAR1(form=~1|TRANSMITTER/YEAR) ,data=data, method="ML")


summary(mod2)


#drop SEX
mod2.a<-update(mod2,.~.-SEX)
#check significance of SEX
anova(mod2,mod2.a)# significant


#drop DAY
mod2.b<-update(mod2,.~.-cDAY)
#check significance of DAY
anova(mod2,mod2.b)# significant!

#drop DAY^2
mod2.c<-update(mod2,.~.-I(cDAY^2))
#check significance of DAY^2
anova(mod2,mod2.c)# insignificant


#drop YEAR
mod2.d<-update(mod2,.~.-YEAR)
#Check significance of YEAR
anova(mod2,mod2.d)# significant difference


#best model

mod.pi<-mod2.c
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
boxplot(log(ACCELMEAN)~SEX,data.ac)
xyplot(log(ACCELMEAN) ~ DAY|TRANSMITTER, type= c("p", "r"), data=data)


#test if random effect better than linear regression
#linear model
gls1<-gls(log(ACCELMEAN)~SEX+cDAY+I(cDAY^2)+YEAR,data=data.ac)
#mixed effects model
mod1<-lme(log(ACCELMEAN)~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER,data=data.ac, method="REML")

#check variance inflation
vif(mod1)

anova(gls1,mod1)#random effect better model, 

#check for autocorrelation
acf(residuals(mod1, type="normalized"))#looks like correlation

#add AR1 correlation structure to model
mod2<-lme(log(ACCELMEAN)~SEX+cDAY+I(cDAY^2)+YEAR, random=~1|TRANSMITTER, correlation=corAR1(form=~1|TRANSMITTER) ,data=data.ac, method="REML")

#the AR1 model is a better fit then the model without the correlation structure
anova(mod1,mod2)

#re-check autocorrelation
acf(residuals(mod2, type="normalized"))#looks pretty good

mod2<-lme(log(ACCELMEAN)~SEX+DAY+I(DAY^2), random=~1|TRANSMITTER, correlation=corAR1(form=~1|TRANSMITTER) ,data=data.ac, method="ML")

#test significance of intereaction
mod2.a<-update(mod2,.~.-SEX)# drop interaction
anova(mod2,mod2.a)#not significant remove

mod2.a<-update(mod2,.~.-DAY)# drop interaction
anova(mod2,mod2.a)#not significant remove

mod2.a<-update(mod2,.~.-I(DAY^2))# drop interaction
anova(mod2,mod2.a)#not significant remove

#best model

mod.ac<-lme(log(ACCELMEAN)~SEX, random=~1|TRANSMITTER, correlation=corAR1(form=~1|TRANSMITTER) ,data=data.ac, method="REML")

r.squaredGLMM(mod2)

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
Anova(mod.ac)


