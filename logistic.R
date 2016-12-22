rm(list=ls())
install.packages(c("elrm","leaps","ResourceSelection","Hmisc"
                   ,"aod","superdiag","mcmcplots","arm"), dependencies 
                 = TRUE, repos = "http://cran.us.r-project.org")
library(R2jags)
library(rjags)
library(runjags)
library(MCMCpack)
library(lattice)
require(elrm)
library(leaps)
library(ResourceSelection)
library(Hmisc)
library(aod)
library(superdiag)
library(mcmcplots)
library(ggmcmc)
library(ggplot2)
library(arm)

dataHDLC=read.csv("~/Documents/memphisclassesbooks/DR. Eugene/dataHDLC.csv")
dataAPOA1=read.csv("~/Documents/memphisclassesbooks/DR. Eugene/dataAPOA1.csv")

fit1=glm(ASCVD~AGE+APO.A1+BMI,family=binomial,data=dataAPOA1)
fit2=glm(ASCVD~AGE+APO.A1,family=binomial,data=dataAPOA1)

summary(fit1)
summary(fit2)

exp(fit1$coef)

fit3=glm(ASCVD~AGE+HDL.C+BMI,family=binomial,data=dataHDLC)
fit4=glm(ASCVD~AGE+HDL.C,family=binomial,data=dataHDLC)

summary(fit3)
summary(fit4)

############################ CIs using profiled log-likelihood  ################

confint(fit2)

##########################################Hosmer and Lemeshow Goodness-of-Fit Test#######################

h2 <- hoslem.test(fit2$y, fitted(fit2), g=8)
h2
#########################odds ratios and 95% CI########################

exp(cbind(OR = coef(fit2), confint(fit2)))


########################## Predict fracture rate for 50 year old########################

pred.age50 <- exp(alpha + b.age*50)/(1+ exp(alpha + b.age*50))




sim.dat <- dataAPOA1   



  bayes.mod <- function(){
    for(i in 1:N){
      
      logit(p[i]) <- alpha + b.age*AGE[i] + b.apoa1*APO.A1[i]+b.bmi*BMI[i]
      ASCVD[i] ~ dbern(p[i])
    }
    alpha ~ dnorm(0.0,1.0E-4) # Prior for intercept
    b.age ~ dnorm(0.0,1.0E-4) # Prior for slope of age
    b.apoa1 ~ dnorm(0.0,1.0E-4) # Prior for slope of APO.A1
    b.bmi~ dnorm(0.0,1.0E-4) # Prior for slope of bmi
  }
  


ASCVD<- dataAPOA1$ASCVD
AGE<-dataAPOA1$AGE
APOA1<-dataAPOA1$APO.A1
BMI <- dataAPOA1$BMI
N <- length(dataAPOA1$AGE)

#sim.dat.jags <- list("ASCVD", "AGE", "APOA1","BMI", "N")
sim.dat.jags <- as.list(dataAPOA1)

sim.dat.jags$N <- nrow(sim.dat)
bayes.mod.params <- c("alpha", "b.age", "b.apoa1","b.bmi")

#bayes.mod.inits <- function(){
#  list("alpha" = 1, "b.age" = 1, "b.apoa1" = 1,"b.bmi"= 1)
#}
bayes.mod.inits <- function(){
  list("alpha" = fit1$coef[1], "b.age" = fit1$coef[2], "b.apoa1" = fit1$coef[3],"b.bmi"= fit1$coef[4])
}





bayes.mod.fit=jags(data = sim.dat.jags, parameters.to.save=bayes.mod.params,
                   inits = bayes.mod.inits, model.file =   bayes.mod)


mod.fit.upd <-update(bayes.mod.fit, n.iter=1000)
bayes.mod.fit.upd <-autojags(bayes.mod.fit)
print(bayes.mod.fit)

par(mfrow=c(2,2))
plot(bayes.mod.fit)
traceplot(bayes.mod.fit)


bayes.mod.fit.mcmc <- as.mcmc(bayes.mod.fit)
summary(bayes.mod.fit.mcmc)

xyplot(bayes.mod.fit.mcmc)


densityplot(bayes.mod.fit.mcmc)


densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

plot(bayes.mod.fit.mcmc)

autocorr.plot(bayes.mod.fit.mcmc)

geweke.diag(bayes.mod.fit.mcmc)

heidel.diag(bayes.mod.fit.mcmc)

############################################
############### ggplots################

install.packages("superdiag", dependencies = TRUE, repos = "http://cran.us.r-project.org")

#library(superdiag)

superdiag(bayes.mod.fit.mcmc, burnin = 100)

install.packages("mcmcplots", dependencies = TRUE, repos = "http://cran.us.r-project.org")

#library(mcmcplots)

denplot(bayes.mod.fit.mcmc)

mcmcplot(bayes.mod.fit.mcmc)

caterplot(bayes.mod.fit.mcmc)


caterplot(bayes.mod.fit.mcmc, parms = c("alpha", "beta1", "beta2"),
          labels = c("alpha", "beta1", "beta2"))



install.packages("ggmcmc", dependencies = TRUE, repos = "http://cran.us.r-project.org")


bayes.mod.fit.gg <- ggs(bayes.mod.fit.mcmc)
ggs_density(bayes.mod.fit.gg)

denplot(bayes.mod.fit.mcmc, parms = c("alpha", "beta1", "beta2"))
traplot(bayes.mod.fit.mcmc, parms = c("alpha", "beta1", "beta2"))


logmcmc = MCMClogit(dataAPOA1$ASCVD~dataAPOA1$AGE+dataAPOA1$APO.A1+dataAPOA1$BMI, burnin=1000,
                    mcmc=21000)

summary(logmcmc)

plot(logmcmc)



logbayes = bayesglm(ASCVD~AGE+APO.A1+BMI,data=dataAPOA1, family=gaussian,prior.mean=0,
                    prior.scale=Inf,prior.df=Inf)

summary(logbayes)

coef=c(fit1$coef[2],fit1$coef[3],fit1$coef[4])
logbayes = bayesglm(ASCVD~AGE+APO.A1+BMI,data=dataAPOA1, family=binomial,prior.mean=coef,
                    prior.scale=Inf,prior.df=Inf)

summary(logbayes)

simulation=coef(sim(logbayes))
par(mfrow=c(2,2))
plot(density(simulation[,1]),ylab="density",xlab="posterior density of Intercept",main="")
plot(density(simulation[,2]),ylab="density",xlab="posterior density of AGE",main="")
plot(density(simulation[,3]),ylab="density",xlab="posterior density of APO.A1",main="")
plot(density(simulation[,4]),ylab="density",xlab="posterior density of BMI",main="")

quantile(simulation[,1],c(.25,.975))
quantile(simulation[,2],c(.25,.975))
quantile(simulation[,3],c(.25,.975))
quantile(simulation[,4],c(.25,.975))