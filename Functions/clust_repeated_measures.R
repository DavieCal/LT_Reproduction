
fit1<-glm(Freq~Var1*Var2, data=clust,family="poisson")
fit2<-glm(Freq~Var1+Var1:Var2, data=clust,family="poisson")
summary(fit1)
summary(fit2)
clust<-rbind(m0,m1,f0,f1)
require(car)
fit1<-Anova(lm(log(Freq)~Var1*Var2,data=clust, contrasts=list(Var1=contr.sum, Var2=contr.sum)), type=3)

fit<-lm(log(Freq)~Var1*Var2,data=clust,contrasts=list(Var1=contr.sum, Var2=contr.sum))
drop1(fit,.~.,test="F")

summary(fit)

summary(aov(log(Freq)~Var1+Var2,data=clust))
TukeyHSD(fit)

require(nlme)
require(multcomp)
require(lsmeans)
require(MuMIn)
mod1<-lme(log10(Freq)~Var1+Var2+Var1:Var2, random=~1|Var3 ,data=clust)
anova(mod1, type="marginal")
mod2<-update(mod1,.~.-Var1:Var2)
anova(mod2, type="marginal")
mod3<-update(mod2,.~.-Var2)
anova(mod3, type="marginal")
summary(mod3)

tuk<-lsmeans(mod3, pairwise~Var1)
r.squaredGLMM(mod3)

summary(x)
anova(mod1,type="marginal")
lmerTest(mod1)
