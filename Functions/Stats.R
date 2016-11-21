#format scientific notation

library(lme4)

data<-daily[daily$MONTH==9,]

#necessary to do likelihood ratio test
fit.full<-lmer(sqrt(STEPLENGTH)~as.factor(SEX)+as.factor(YEAR)+(1|TRANSMITTER),data=data, REML=FALSE)
fit.2<-lmer(sqrt(STEPLENGTH)~as.factor(SEX)+(1|TRANSMITTER),data=data, REML=FALSE)
fit.1<-lmer(sqrt(STEPLENGTH)~as.factor(YEAR)+(1|TRANSMITTER),data=data, REML=FALSE)
fit.null<-lmer(sqrt(STEPLENGTH)~1+(1|TRANSMITTER),data=data, REML=FALSE)

fit.full.slope<-lmer(sqrt(STEPLENGTH)~as.factor(SEX)+as.factor(YEAR)+(1+YEAR|TRANSMITTER),data=data, REML=FALSE)
#look at residuals
plot(fit.full.slope)
#normality of residuals
hist(residuals(fit.full.slope))
fit<-lm(sqrt(STEPLENGTH)~as.factor(SEX)+as.factor(YEAR),data=data)

boxplot(STEPLENGTH~YEAR*SEX,data)
head(dfbeta(fit))

