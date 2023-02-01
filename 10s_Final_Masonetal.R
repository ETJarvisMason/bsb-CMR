# Mason et al.
# Capture-mark-reencounter (CMR) model for modeling L. Bellquist 2010s CMR data
# 2022

# Notes on L. Bellquist BSB CMR data--------------------------------------------

# Tagging occurred largely in San Diego area, not throughout so. California
# Data for each fish includes whether fish were kept or released alive (w/wo tag)
# Some fish were tagged and released in Marine Protected Areas (MPAs), HOWEVER
#   Only fish released outside of MPAs are included
# More consistent monthly sampling
#   Monthly occasions from Oct 2012 through Jan 2015
#   Tag loss rate adjusted to monthly rate in model
# There is only 1 recapture/recovery with associated length data
# Release lengths were reported in standard length (SL)
#   SL converted to total length using conversion in Love et al. 1997
# Annual phi estimated in derived quantities section
# The 2 inch increase in the Minimum Size Limit (MSL) took effect 
#   March 2013, thus the code for assigning the size classification was
#   updated to use the old MSL for occs 1-5 and the new MSL for occs > 5

# load libraries
library(tidyverse)
library(jagsUI)

# Define model for 2010s CMR data--------------------------------------------
sink("tvr_size10s.jags")
cat("
  model{

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # priors
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # annual tag retention rate
    r ~ dbeta(140.2,26.8) # > rb_parms:  $alpha [1] 140.1815 $beta [1] 26.81656
    r.mo <- r^(1/12) # monthly rate
    
    # assign resighting and kill prob according to size class and season
    for (t in 1:(nT-1)){
     
     phi[1,t] <- phi.grp[1,seas[t+1]] 
     phi[2,t] <- phi.grp[2,seas[t+1]] 
     
     R[1,t] <- R.grp[1,seas[t+1]]          # sublegal by season (spawn/non-spawn)
     R[2,t] <- R.grp[2,seas[t+1]]          # legal by season (spawn/non-spawn)
     
     kappa[1,t] <- kappa.grp[1,seas[t]]
     kappa[2,t] <- kappa.grp[2,seas[t]]
     }#t
     
    phi.grp[1,1] ~ dbeta(1,1)
    phi.grp[1,2] ~ dbeta(1,1)
     
    phi.grp[2,1] ~ dbeta(1,1)
    phi.grp[2,2] ~ dbeta(1,1)
     
    R.grp[1,1] ~ dbeta(1,1)
    R.grp[1,2] ~ dbeta(1,1)
    
    R.grp[2,1] ~ dbeta(1,1)
    R.grp[2,2] ~ dbeta(1,1)
    
    kappa.grp[1,1] ~ dbeta(1,1)
    kappa.grp[1,2] ~ dbeta(1,1)
    
    kappa.grp[2,1] ~ dbeta(1,1)
    kappa.grp[2,2] ~ dbeta(1,1)
    
    for (t in 1:(nT-1)){          # for every occasion
      p[t] <- p.occ[occ_type[t]]  # select p prob according to survey/no survey
    }#t

    p.occ[1] ~ dbeta(1,1)
    p.occ[2] <- 0    
    
    # beta distributions defining hypothetical tag reporting rates of 25, 50, and 75%
    theta[1] ~ dbeta(75,25)
    theta[2] ~ dbeta(50,50)
    theta[3] ~ dbeta(25,75)

    # for t1 = age 24 mo, t3 = age 192 mo
    L1 ~ dnorm(191,100)
    L2 ~ dnorm(391,100)
    L3 ~ dnorm(487,100)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Growth
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 1:nI){
      
      s[i,f[i]] <- s.dat[i,f[i]]           # size at tagging assigned from s.dat
      
      for (t in (f[i] + 1):nT){
        
        s[i,t] <- s[i,t-1] +
                  (-2 * (log((L3 - L2) / (L2 - L1))) / (t3 - t1)) *
                  ((L1 + ((L3 - L1) / (1 - ((L3 - L2) / (L2 - L1))^2)))
                  - s[i,t-1])
        
      s.dat[i,t] ~ dnorm(s[i,t], 100)    # where fish size data are filled in...drawn from a normal distribution of size s at time t, where s is being estimated based on growth
      }#t
      
      
      # for fish growing through occ 5 (through Feb 2013), use MSL 305 to assign class, otherwise use MSL 356
      for (t in f[i]:nT){
      class[i,t] <- ifelse(t < 6, ifelse(s[i,t] >= 305, 2, 1),
                                  ifelse(s[i,t] >= 356, 2, 1))  # assign class based on length at time t
      }#t
    }#i

    #--------------------------------------
    # Latent states (s):
    # 1) alive, with tag
    # 2) caught & kept
    # 3) unavailable
    #--------------------------------------


    for (i in 1:nI){
      for (t in f[i]:(nT-1)){
        ps[1,1,i,t] <- phi[class[i,t],t] * r.mo
        ps[1,2,i,t] <- kappa[class[i,t],t] * r.mo
        ps[1,3,i,t] <- (1 - phi[class[i,t],t] - kappa[class[i,t],t]) * r.mo + (1 - r.mo)

        ps[2,1,i,t] <- 0
        ps[2,2,i,t] <- 0
        ps[2,3,i,t] <- 1

        ps[3,1,i,t] <- 0
        ps[3,2,i,t] <- 0
        ps[3,3,i,t] <- 1

    #--------------------------------------
    # Tagged fish observations (o):
    # 1) recap by biologist and angler
    # 2) recap by biologist only
    # 3) recap by angler only
    # 4) caught and kept by angler
    # 5) not caught or reported
    #--------------------------------------
    


        po[1,1,i,t] <- p[t] * R[class[i,t],t]
        po[1,2,i,t] <- p[t] * (1 - R[class[i,t],t])
        po[1,3,i,t] <- (1 - p[t]) * R[class[i,t],t]
        po[1,4,i,t] <- 0
        po[1,5,i,t] <- (1 - p[t]) * (1 - R[class[i,t],t])

        po[2,1,i,t] <- 0
        po[2,2,i,t] <- 0
        po[2,3,i,t] <- 0
        po[2,4,i,t] <- 1
        po[2,5,i,t] <- 0

        po[3,1,i,t] <- 0
        po[3,2,i,t] <- 0
        po[3,3,i,t] <- 0
        po[3,4,i,t] <- 0
        po[3,5,i,t] <- 1

     }#t
   }#i

  # marginalized capture recapture likelihood
  for (i in 1:nI){
    zeta[i,f[i],1] <- 1
    zeta[i,f[i],2] <- 0
    zeta[i,f[i],3] <- 0

    for(t in f[i]:(nT-1)){
      zeta[i,(t+1),1] <- inprod(zeta[i,t,], ps[,1,i,t]) * po[1,y[i,t+1],i,t]
      zeta[i,(t+1),2] <- inprod(zeta[i,t,], ps[,2,i,t]) * po[2,y[i,t+1],i,t]
      zeta[i,(t+1),3] <- inprod(zeta[i,t,], ps[,3,i,t]) * po[3,y[i,t+1],i,t]

      }#t
      
      lik[i] <- sum(zeta[i, nT, ])
      one[i] ~ dbern(lik[i])
  }#i

# derived quantities

for (i in 1:2){
ann_phi[i] <- phi.grp[i,1]^2 * phi.grp[i,2]^9
}#i

month.k[1,1] <- kappa.grp[1,1] # sublegal summer kappa month 1
month.k[2,1] <- kappa.grp[2,1] # legal summer kappa month 1

for (i in 1:2){
for (t in 2:2){
  month.k[i,t] <- (phi.grp[i,1]^(t-1)) * kappa.grp[i,1]
}#t
}#i

for (i in 1:2){
for (t in 3:11){
  month.k[i,t] <- (phi.grp[i,1]^2) * (phi.grp[i,2]^(t-3)) * kappa.grp[i,2]
  }#t
}#i
for (i in 1:2){
ann.kappa[i] <- sum(month.k[i,])
}#i

# Pop estimates assuming 100% tag reporting
for (t in 1:nH){

    N.legal[t] <- (H.legal[t] / ann.kappa[2]) / 1000 # report in millions of fish
    N.sublegal[t] <- (H.sublegal[t] / ann.kappa[1]) / 1000 # report in millions of fish
    N.tot[t] <- N.legal[t] + N.sublegal[t] # in millions of fish
}#t

# conditional population estimates (according to hypothetical tag reporting (theta) rates)
for (t in 1:nH){
    for (j in 1:3){
    N.legal.theta[t,j] <- (H.legal[t] / ann.kappa[2] * theta[j]) / 1000 # report in millions of fish
    N.sublegal.theta[t,j] <- (H.sublegal[t] / ann.kappa[1] * theta[j]) / 1000 # report in millions of fish
    N.tot.theta[t,j] <- N.legal.theta[t,j] + N.sublegal.theta[t,j] # in millions of fish
}#j
}#t


}",fill = TRUE)
sink()

# Load L. Bellquist BSB CMR data (2012-2015)-------------------------------------------------------------------------

# load encounter history and length data
load("ch10sAnn.Rdata")
load("ch10sLen.Rdata")

df <- arrange(dat.final,fish2) # ch obs matrix

y2 <- data.matrix(df[,2:29]) # get rid of extra columns 

# Compute date of first capture
get.first <- function(x) min(which(x!=5))
f <- apply(y2, 1, get.first)

# starts October 2012,ends January 2015 ( spawning season is June, July, August )
# there are 28 occs (months)
seas <- c(rep(2,8),rep(1,3),rep(2,9),rep(1,3),rep(2,6)) # season of month (1=spawn,2=non-spawn)

s.dat2 <- data.matrix(l[,-c(1:2)]) # remove the 2 fish id columns

hist(s.dat2)

one2 <- rep(1, nrow(y2))

# create survey vector, call it occ_type (1 = survey, 0 = no survey)
nT <- ncol(y2)
occ_type <- as.vector(matrix(2,nrow=nT-1))
sampoccs <- unique(f)      # just survey occs
occ_type[sampoccs] <- 1
occ_type[is.na(occ_type)] <- 2

# Set up data ----------------------------------------------------

# Bring in total harvest estimates
load("BSB_est_harvest.Rdata")

# harvest data represent thousands of fish
H.legal <- H.legal.10s
H.sublegal <- H.sublegal.10s

nH = length(H.legal)

# ages for Francis parameterization of Von Bert growth
t1 <- 24    # fixed, age class 2
t3 <- 192   # fixed, age class 16

# Data defined for model
jags.data <- list(y = y2, nI = nrow(y2), nT = ncol(y2), f = f, one = one2, s.dat = s.dat2,
                  t1 = t1, t3 = t3, seas = seas, occ_type = occ_type,
                  H.legal = H.legal, H.sublegal = H.sublegal, nH = nH)

# Parameters monitored
parameters <- c('phi.grp','kappa.grp','R.grp','p.occ',
                'ann_phi','ann.kappa', 'theta',
                'L1','L2','L3',
                'N.tot','N.tot.theta')

# MCMC specification
nc <- 3
nt <- 5
ni <- 20000
nb <- 5000

# Run CMR Model in JAGS ------------------- 

set.seed(2022)

jags.10s <- jags(jags.data, inits = NULL, parameters, "tvr_size10s.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb,
                parallel=T, n.cores = nc)

moddat <- jags.10s

library(MCMCvis)

params <- parameters

MCMCtrace(moddat,
          pdf = TRUE,
          params = params,
          Rhat = TRUE,
          n.eff = TRUE,
          post_zm = TRUE)

save(jags.10s,file="10s_final.Rdata")

 
