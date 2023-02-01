# Mason et al.
# Age and Growth Model Fits (Traditional Von Bertalanffy and Francis parameterization)
# 2022

# load libraries

library(FSA)
library(nlstools)
library(readr)

# load data (Barred Sand Bass age at length data from Walker et al. 2020)
bsb <- read.csv("bsbage.csv")
str(bsb)

# create column of age in months (for CMR modeling of 2010s tagging data)
bsb.mo <- bsb %>% 
  mutate(age.mo = age*12)
View(bsb.mo)

hist(bsb$age)
hist(bsb.mo$age.mo)


# Annual Growth -----------------------------------------------------------

# traditional Von Bert parameterization
svTypical <- vbStarts(tl~age,data=bsb)
unlist(svTypical) 

vbTypical <- tl~Linf*(1-exp(-K*(age-t0)))
fitTypical <- nls(vbTypical,data=bsb,start=svTypical)

fitPlot(fitTypical,xlab="Age",ylab="Total Length (mm)",main="")
overview(fitTypical)

# note high correlation among traditional Von Vertalanffy growth parameters
bootTypical <- nlsBoot(fitTypical,niter=1000)
confint(bootTypical,plot=TRUE)
plot(bootTypical)


# Francis parameterization of the Von Bertalanffy Growth Function
svFrancis <- vbStarts(tl~age,data=bsb,type="Francis",tFrancis=c(2,9))
unlist(svFrancis)

table(bsb$age)
ages <- c(3,16)

vbFrancis <- vbFuns("Francis")
fitFrancis <- nls(tl~vbFrancis(age,L1,L2,L3,t1=ages[1],t3=ages[2]),
                  data=bsb,start=svFrancis)
overview(fitFrancis)

# solves correlation problem
bootFrancis <- nlsBoot(fitFrancis,niter=1000)
confint(bootFrancis,plot=TRUE)
plot(bootFrancis)

# Annual Growth -----------------------------------------------------------

# traditional Von Bert parameterization
svTypical <- vbStarts(tl~age.mo,data=bsb.mo)
unlist(svTypical) 

vbTypical <- tl~Linf*(1-exp(-K*(age.mo-t0)))
fitTypical <- nls(vbTypical,data=bsb.mo,start=svTypical)

fitPlot(fitTypical,xlab="Age (months)",ylab="Total Length (mm)",main="")
overview(fitTypical)

bootTypical <- nlsBoot(fitTypical,niter=1000)
confint(bootTypical,plot=TRUE)
plot(bootTypical)

# Francis parameterization of the Von Bertalanffy Growth Function
svFrancis <- vbStarts(tl~age.mo,data=bsb.mo,type="Francis",tFrancis=c(24,192))
unlist(svFrancis)

table(bsb.mo$age.mo)
ages.mos <- c(24,192)

vbFrancis <- vbFuns("Francis")
fitFrancis <- nls(tl~vbFrancis(age.mo,L1,L2,L3,t1=ages.mos[1],t3=ages.mos[2]),
                  data=bsb.mo,start=svFrancis)
overview(fitFrancis)

bootFrancis <- nlsBoot(fitFrancis,niter=1000)
confint(bootFrancis,plot=TRUE)
plot(bootFrancis)


