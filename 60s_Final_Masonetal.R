# Mason et al.
# Capture-mark-reencounter (CMR) model for modeling the 1960s CDFW CMR data
# 2022

# load libraries
library(tidyverse)
library(jagsUI)

# Define model for 1960s CDFW CMR data--------------------------------------------

sink("tvr_size60s.jags")
cat("
    
model{    
    phi[1] ~ dbeta(1,1)
    phi[2] ~ dbeta(1,1)
    
    R[1] ~ dbeta(1,1)
    R[2] ~ dbeta(1,1)
    
    kappa[1] ~ dbeta(1,1)
    kappa[2] ~ dbeta(1,1)

    p ~ dbeta(1,4)
    
    # beta distributions defining hypothetical tag reporting rates of 25, 50, and 75%
    theta[1] ~ dbeta(75,25)
    theta[2] ~ dbeta(50,50)
    theta[3] ~ dbeta(25,75)
    
    # for t1 = age 3, t3 = age 16
    L1 ~ dnorm(236,10)
    L2 ~ dnorm(403,10)
    L3 ~ dnorm(495,10)

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
      for (t in f[i]:nT){
        class[i,t] <- ifelse(s[i,t] >= 305, 2, 1) # assign class based on length at time t
      }#t
    }#i  
  
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # tag retention as a function of years at liberty (age of tag)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (j in 1:nT) {
      tr[j] ~ dbeta(trbeta_a[j], trbeta_b[j])
    }#j
    
    #--------------------------------------
    # Latent states (s):
    # 1) alive, with tag
    # 2) caught & kept
    # 3) unavailable
    #--------------------------------------

    for (i in 1:nI){
      for (t in f[i]:(nT-1)){
        
        ps[1,1,i,t] <- phi[class[i,t]] * tr[t-f[i] + 1]
        ps[1,2,i,t] <- kappa[class[i,t]] * tr[t-f[i] + 1]
        ps[1,3,i,t] <- (1 - phi[class[i,t]] - kappa[class[i,t]]) * tr[t-f[i] + 1] + (1 - tr[t-f[i] + 1])
    
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
        po[1,1,i,t] <- p * R[class[i,t]]
        po[1,2,i,t] <- p * (1 - R[class[i,t]])
        po[1,3,i,t] <- (1 - p) * R[class[i,t]]
        po[1,4,i,t] <- 0
        po[1,5,i,t] <- (1 - p) * (1 - R[class[i,t]])    
    
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
  
  # estimated population size of legal/sublegal fish
  # for the 1960s, our best estimate of total harvest is from 1964

  for (t in 1:nH){

    # assuming 100% tag reporting
    N.legal[t] <- (H.legal[t] / kappa[2]) / 1000 # report in millions of fish
    N.sublegal[t] <- (H.sublegal[t] / kappa[1]) / 1000 # report in millions of fish
    N.tot[t] <- N.legal[t] + N.sublegal[t] # in millions of fish
    
    # conditional population estimates (according to hypothetical tag reporting (theta) rates)
    for (j in 1:3){
    
    N.legal.theta[t,j] <- H.legal[t] / (kappa[2] / theta[j]) / 1000 # report in millions of fish
    N.sublegal.theta[t,j] <- H.sublegal[t] / (kappa[1] / theta[j]) / 1000 # report in millions of fish
    N.tot.theta[t,j] <- N.legal.theta[t,j] + N.sublegal.theta[t,j] # in millions of fish
    }#j  
  }#t

}",fill = TRUE)
sink()

# Load CDFW BSB CMR data (1962-1970)-------------------------------------------------------------------------

# load encounter history and length data
load("ch60sAnn2.Rdata")
load("ch60sLen2.Rdata")

df <- dat.final # ch obs matrix

y2 <- data.matrix(df[,2:10]) # get rid of extra columns 

# Compute date of first capture
get.first <- function(x) min(which(x!=5))
f <- apply(y2, 1, get.first)

# length matrix
s.dat2 <- data.matrix(l[,-1]) # remove fish id column

hist(s.dat2)

one2 <- rep(1, nrow(y2))


# Set up data ----------------------------------------------------

# Bring in total harvest estimates
load("BSB_est_harvest.Rdata")

# harvest data represent thousands of fish
H.legal <- H.legal.64
H.sublegal <- H.sublegal.64

nH = length(H.legal)

# year-specific tag retention beta distribution params for defining tag retention priors
load("trbeta_parms.Rdata")

# ages for Francis parameterization of Von Bert growth
t1 = 3
t3 = 16


# Data defined for model
jags.data <- list(y = y2, nI = nrow(y2), nT = ncol(y2), f = f, one = one2, s.dat = s.dat2, t1 = t1, t3 = t3, 
                  trbeta_a = trbeta_a, trbeta_b = trbeta_b,
                  H.legal = H.legal, H.sublegal = H.sublegal, nH = nH)

# Parameters to monitor
parameters <- c('phi','kappa','R','p','theta',
                'L1','L2','L3',
                'N.tot','N.tot.theta')

# MCMC specification
nc <- 3
nt <- 5
ni <- 20000
nb <- 5000

# Run CMR Model in JAGS ------------------- 

set.seed(2022)

jags.60s <- jags(jags.data, NULL, parameters, "tvr_size60s.jags", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb,
                  parallel=T, n.cores = nc)

moddat <- jags.60s

library(MCMCvis)

print(moddat,digits = 3)

params <- parameters

MCMCsummary(moddat,
            params = params)

MCMCtrace(moddat,
          pdf = TRUE,
          params = params,
          Rhat = TRUE,
          n.eff = TRUE,
          post_zm = TRUE)

save(jags.60s,file="60s_final.Rdata")