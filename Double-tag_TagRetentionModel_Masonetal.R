# Mason et al.
# Tag retention model for double-tagging data (2 identical tags)
# with individual recovery times
# 2022

# Background and data simulations
# Explore relationship between the ratio of fish with 2 vs 1 tag--------

tagpop<-1000 # num of tagged fish
fish<-matrix(data=NA, nrow = tagpop,ncol=3) # matrix for storing sims
ratioStor<-matrix(data=NA, nrow = 50,ncol=2) # matrix for storing ratio of tags
taglossrate<-0 # init for the for loop below

# over 50 time steps, with a 0.015 rate of tag loss at each step, calculate
for (i in 1:50){ 
  taglossrate <- taglossrate + .015
  
  fish[,1]<-rbinom(tagpop, 1, (1-taglossrate)) # prob of keeping first tag
  fish[,2]<-rbinom(tagpop, 1, (1-taglossrate)) # prob of keeping second tag
  fish[,3]<-fish[,1]+fish[,2] # for each fish, calculate the number of tags remaining after a single time step
  
  ratioStor[i,1]<-length(which(fish[,3]==1)) / length(which(fish[,3]==2)) # What is the ratio of single to double tagged recaps?
  ratioStor[i,2]<- taglossrate
}

plot(ratioStor[,2],ratioStor[,1]) # plot ratio of 1:2 tags retained in recaps against tag loss rate

# Note important relationship between the odds of losing a tag and the ratio of 1:2 tags retained
odds<-ratioStor[,2] / (1-ratioStor[,2]) *2 
# odds of losing a tag = the ratio of the probability of losing a tag to the prob. of retaining a tag (x2 since we have two tags);
# explains the ratio of 1:2 tags retained in recaps!!
plot(odds,ratioStor[,1]) # note direct proportionality

# Next steps:
# We have recaps by time, so we need to:
# A) Calculate the time-specific prob. of tag loss, and
# B) Calculate the prob. that a double tagged fish has retained one vs. two tags (see above)


# Simulate data and define model --------------------------------------------

# sample recapture times over a 2 year window:
capT<-runif(tagpop, min=0.000001, max=3)
# create new matrix to store simulated fish with tags
fishSim<-matrix(data=NA, nrow = tagpop,ncol=3)

# For each of those fish, calculate the prob that a fish loses 1 or 2 tags based on recap time, as a function of age of tag (time at liberty)
# note: log is the natural log in R

# Assume yearly tag retention rate of 80%
trr<-log(0.8) # continuous rate of loss in log space

# Assume exponential decay for tag loss
qt<-2.5       # quadratic term for capturing increase in prob. of loss over time

# Incorporate immediate tag loss
itl<-0.98     # probability of retaining tag immediately after release

# Simulate tag loss of first tag for each individual based on recap time
fishSim[,1]<-rbinom(tagpop,1,itl*exp(trr*capT^qt)) 
# Next simulate tag loss of 2nd tag for each individual based on recap time
fishSim[,2]<-rbinom(tagpop,1,itl*exp(trr*capT^qt))

# How many tags does each fish have?
fishSim[,3]<-fishSim[,1]+fishSim[,2]

# Remove all fish that lost both tags from our data set
OnlyRecap<-fishSim[which(fishSim[,3]>0),3]
RecapT<-capT[which(fishSim[,3]>0)] # keep only the recap times of fish with at least one tag remaining

# Finally, recode our data for binomial likelihood such that:
# a) a fish with one tag is coded as 0
recode_d<-OnlyRecap*0
# b) a fish with two tags is coded as 1
recode_d[which(OnlyRecap==2)]<-1
N<-length(recode_d) # length of our recapture data

# fit a model in JAGS to make sure we can recover tag loss...

# libraries for running JAGS
library(runjags)
library(rjags)
library(R2jags)

# put data into JAGS readble format
dat <- list(recode_d= recode_d, N=N, RecapT=RecapT)

#Tag retention model (tag loss as a function of age of tag)

tag.model <- "model {

# type I loss (immediate); 1 - alpha
alpha ~ dbeta(1,1)   # probability of retaining a tag immediately after release

# type II tag loss (annual chronic (continuous) rate of loss)
beta ~ dnorm(0,0.01) # also equal to the natural log of the annual discrete rate of tag retention 
RR <- exp(-beta)     # annual discrete rate of tag retention 

# quadratic term to account for increase in tag loss prob. with age of tag 
lambda ~ dlnorm(0,1)


for (i in 1:N){

# for each fish,
# calculate the probability it still has a tag at time t

Q[i] <- alpha * exp(-(beta * RecapT[i]^lambda))

# fish-specific probability of losing a tag
pf[i] <- 1 - Q[i]

# From exercise above, we know that the ratio of 1:2 tags among recaptures is equal to the odds of losing a tag:
# P(1) / (1-P(1)) = R / (1- R) *2
# ... where P(1) is the prob of a fish having one tag (vs. 2 tags -- {1-P(1)})
# ... and where R is the tag loss probability

# to convert from odds to proportion --> odds / (1+odds)
# Thus, P(1) = [R / (1-R) *2] / 1+(R / (1-R) *2)]

# Calculate probability of losing one tag
P1[i]<- (pf[i] / (1-pf[i]) *2) / (1+(pf[i] / (1-pf[i]) *2)) 
# Calculate probability of retaining both tags (not losing a tag)
p2[i]<- 1-P1[i]

# Note also the following is true:
# lost both tags
# p0[i] <- (1 - Q[i]) * (1 - Q[i])
# kept tag 1 and lost tag 2 or lost tag 1 and kept tag 2
# p1[i] <- (1 - Q[i]) * Q[i] + Q[i] * (1 - Q[i])    
# kept both tags
# p2[i] <- Q[i] * Q[i]

} #end loop

# Likelihood
for (i in 1:N){
recode_d[i] ~ dbin(p2[i],1)
}#i

}"


# Test model --------------------------------------------------------------


tag.out <- run.jags(tag.model, monitor=c("deviance", 'alpha', 'beta','lambda','RR'), 
                    data=dat, burnin = 500, sample = 1000, 
                    n.chains=3)

# Model works...recovers tag loss parameters with uncertainty
# Note that trr = -beta, alpha = itl, lamda = qt  

summary(tag.out)
plot(tag.out, vars = c('alpha', 'beta','lambda','RR'))

# Load and format Kelp Bass double tag data ----------------------------------------------------

library(tidyverse)
library(lubridate)
library(jagsUI)

kbtags <- read.csv("KB_Double_tagging_data_etm.csv",stringsAsFactors = FALSE)
names(kbtags)
str(kbtags)

# data formatting
data <- kbtags %>% 
  dplyr::select(Tag.ID.1,Date.Tagged,R1Date,R2Date,R3Date,R1Tag,R2Tag,R3Tag) %>%
  mutate(tagdate = mdy(Date.Tagged), rdate1 = mdy(R1Date), rdate2 = mdy(R2Date),
         rdate3 = mdy(R3Date)) %>% 
  dplyr::select(Tag.ID.1,tagdate,rdate1,rdate2,rdate3,R1Tag,R2Tag,R3Tag) %>% 
  mutate(dal1 = rdate1 - tagdate, dal2 = rdate2 - tagdate, dal3 = rdate3 - tagdate)

# first-time recaptures with 2 tags
recap1_2 <- data %>% 
  dplyr::select(R1Tag,dal1)%>% 
  dplyr::filter(R1Tag == 2) %>% 
  dplyr::rename(DAL = dal1,BB = R1Tag)
head(recap1_2)

# second-time recaptures with 2 tags
recap2_2 <-  data %>% 
  dplyr::select(R2Tag,dal2)%>% 
  dplyr::filter(R2Tag == 2)%>% 
  dplyr::rename(DAL = dal2,BB = R2Tag)
head(recap2_2) 

# third-time recaptures with 2 tags
recap3_2 <-  data %>% 
  dplyr::select(R3Tag,dal3)%>% 
  dplyr::filter(R3Tag == 2)%>% 
  dplyr::rename(DAL = dal3,BB = R3Tag)
head(recap3_2) 

# create a single table of fish retaining 2 tags
BB_BB <- rbind(recap1_2,recap2_2,recap3_2) 
BB_BB$Type <- "AA"

# first-time recaptures with only 1 tag remaining
recap1_1 <- data %>% 
  dplyr::select(R1Tag,dal1)%>% 
  dplyr::filter(R1Tag == 1)%>% 
  dplyr::rename(DAL = dal1,BB = R1Tag) %>% 
  drop_na()
head(recap1_1)

# second-time recaptures with only 1 tag remaining
recap2_1 <-  data %>% 
  dplyr::select(R2Tag,dal2)%>% 
  dplyr::filter(R2Tag == 1)%>% 
  dplyr::rename(DAL = dal2,BB = R2Tag)
head(recap2_1) 

# third-time recaptures with only 1 tag remaining
recap3_1 <-  data %>% 
  dplyr::select(R3Tag,dal3)%>% 
  dplyr::filter(R3Tag == 1)%>% 
  dplyr::rename(DAL = dal3,BB = R3Tag)
head(recap3_1) 

# create a single table of fish retaining 1 tag
BB_B <- rbind(recap1_1,recap2_1,recap3_1) 
BB_B$Type <- "A"


# combine single and double tag returns into single dataframe
df <- rbind(BB_BB,BB_B) %>%  
  arrange(DAL) 
df$DAL <- as.numeric(df$DAL)
hist(df$DAL)

# convert data to format needed to run Brice's JAGS code
OnlyRecap <- df$BB
RecapT <- df$DAL/365 # convert to years

RecapT

recode_d<-OnlyRecap*0 # only has 1 tag
recode_d[which(OnlyRecap==2)]<-1 # still has 2 tags
N<-length(recode_d) 


#put data into JAGS readble format
dat <- list(recode_d= recode_d, N=N, RecapT=RecapT)


# Tag retention model (tag loss as a function of age of tag) -------


sink("tag_loss_tvr.jags")
cat(" 
model{

# type I loss (immediate); 1 - alpha

alpha ~ dbeta(1,1)   # probability of retaining a tag immediately after release

# type II tag loss (annual continuous rate of loss)

beta ~ dnorm(0,0.01) # also equal to the natural log of the annual discrete rate of tag retention 
RR <- exp(-beta)     # annual discrete rate of tag retention 

# quadratic term to account for increase in tag loss prob. with age of tag 
lambda ~ dlnorm(0,1)


# for each fish, calculate the probability of retaining a tag to time t

for (i in 1:N){

Q[i] <- alpha * exp(-(beta * RecapT[i]^lambda))

pf[i] <- 1 - Q[i] # probability of losing a tag

# From exercise above, we know that the ratio of 1:2 tags among recaptures is equal to the odds of losing a tag:
# P(1) / (1-P(1)) = R / (1- R) *2
# ... where P(1) is the prob of a fish having one tag (vs. 2 tags -- {1-P(1)})
# ... and where R is the tag loss probability

# to convert from odds to proportion --> odds / (1+odds)
# Thus, P(1) = [R / (1-R) *2] / 1+(R / (1-R) *2)]

# Calculate probability of losing one tag
P1[i]<- (pf[i] / (1-pf[i]) *2) / (1+(pf[i] / (1-pf[i]) *2)) 
# Calculate the probability of losing tags
p2[i]<- 1-P1[i] 

} #end loop

# Likelihood
for (i in 1:N){
recode_d[i] ~ dbin(p2[i],1)
}#i

}",fill = TRUE)
sink()

params <- c("deviance", 'alpha', 'beta','lambda','RR')


# Run tag retention model on double-tag data ------------------------------

set.seed(010203)

tag.out <- jags(model = "tag_loss_tvr.jags", 
                parameters = params,
                data=dat, n.burnin = 5000, n.iter = 10000, 
                n.chains=3)

print(tag.out)
plot(tag.out)


# Plot cumulative and non-cumulative tag retention over time ---------------------------------


samps <- length(tag.out$sims.list$alpha)
max.age <- 11
retain <- matrix(NA, samps, max.age)
lost <- matrix(NA, samps, max.age)


# cumulative proportion of recaptured fish with a tag
for (i in 1:samps){
  for (j in 1:max.age){
    retain[i,j] <- tag.out$sims.list$alpha[i] * exp(-(tag.out$sims.list$beta[i]*j^tag.out$sims.list$lambda[i]))
  }
}

# probability of tag retention, as a function of time at liberty
lost[,1] <- retain[,1]
for (j in 2:max.age){
  lost[,j] <- 1- (retain[,(j-1)] - retain[,j])/(retain[,(j-1)])
}

par(family = 'serif', mfrow = c(1,1))
# plot cumulative tag retention probability 
boxplot(retain, ylim = c(0,1))

cum.prob.tab <- gather(as.data.frame(retain)) %>% 
  rename(TAL = key, prob = value) %>% 
  mutate(yr = case_when(TAL == "V1" ~ 1,
                        TAL == "V2" ~ 2,
                        TAL == "V3" ~ 3,
                        TAL == "V4" ~ 4,
                        TAL == "V5" ~ 5,
                        TAL == "V6" ~ 6,
                        TAL == "V7" ~ 7,
                        TAL == "V8" ~ 8,
                        TAL == "V9" ~ 9,
                        TAL == "V10" ~ 10,
                        TAL == "V11" ~ 11))

cum <- ggplot(cum.prob.tab, aes(x=as.factor(yr), y=prob))+
  geom_boxplot(outlier.shape = NA) + 
  xlab("Time at liberty (years)")+
  ylab("Cumulative proportion")+
  theme_classic()
cum

ret.prob.tab <- gather(as.data.frame(lost)) %>% 
  rename(TAL = key, prob = value) %>% 
  mutate(yr = case_when(TAL == "V1" ~ 1,
                        TAL == "V2" ~ 2,
                        TAL == "V3" ~ 3,
                        TAL == "V4" ~ 4,
                        TAL == "V5" ~ 5,
                        TAL == "V6" ~ 6,
                        TAL == "V7" ~ 7,
                        TAL == "V8" ~ 8,
                        TAL == "V9" ~ 9,
                        TAL == "V10" ~ 10,
                        TAL == "V11" ~ 11))


noncum <- ggplot(ret.prob.tab, aes(x=as.factor(yr), y=prob))+
  geom_boxplot(outlier.shape = NA) + 
  xlab("Time at liberty (years)")+
  ylab("Retention probability")+
  theme_classic()
noncum

# plot year-specific probability (non-cum.) of tag retention
boxplot(lost)

# calculate mean year-specific tag retention probs (non-cum.)
mu <- rep(NA, max.age)
sig <- rep(NA, max.age)

for (j in 1:max.age){
  mu[j] <- mean(lost[,j],na.rm = TRUE)
  sig[j] <- sd(lost[,j],na.rm = TRUE)  
}

# calculate beta distribution shape parameter values
# for use as year-specific tag retention priors in CMR model 
trbeta_a <- mu * 1/(sig * sig)
trbeta_b <- (1-mu) * 1/(sig*sig)

save(trbeta_a,trbeta_b,file = "trbeta_parms.Rdata")

# check the distributional form of each year-specific beta distribution
a <- trbeta_a
b <- trbeta_b

hist(rbeta(10000,a[1],b[1]))
hist(rbeta(10000,a[2],b[2]))
hist(rbeta(10000,a[3],b[3]))
hist(rbeta(10000,a[4],b[4]))
hist(rbeta(10000,a[5],b[5]))
hist(rbeta(10000,a[6],b[6]))
hist(rbeta(10000,a[7],b[7]))
hist(rbeta(10000,a[8],b[8]))
hist(rbeta(10000,a[9],b[9]))
hist(rbeta(10000,a[10],b[10]))
hist(rbeta(10000,a[10],b[11]))
