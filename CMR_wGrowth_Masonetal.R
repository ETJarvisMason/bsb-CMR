# Mason et al.
# Capture-mark-reencounter (CMR) data simulation and model testing
# 2022

# load libraries
library(tidyverse)
library(jagsUI)
library(MCMCvis)


# Simulate CMR data -------------------------------------------------------


set.seed(1236)

nT <- 5                           # occasions
rel <- rep(50,nT-1)               # individuals released per year
f <- rep(1:(nT-1), times = rel)   # first release occasion
nI = length(f)                    # total number of individuals released
limit <- 305                      # 305 mm size (harvest) limit

# size of each individual at release
sd.s <- 3
s <- rnorm(nI, 320, 75)
par(family = 'serif', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1))
hist(s, main = NULL, xlab = 'Length (cm)', cex.lab = 2, breaks = 25)

phi <- 0.75          # survival probability
R <- c(0.1,0.01)     # angler catch-and-release probability
k <- c(0.01,0.15)    # recovery (harvest) probability or kill (k) rate

p <- 0.2             # biologist recapture probability
r <- 0.84            # annual tag retention rate, here it is a constant rate

k.mat <- matrix(NA, nI, nT-1)
R.mat <- matrix(NA, nI, nT-1)
s.mat <- matrix(NA, nI, nT)

Psi <- array(NA, dim = c(3,3,nI,nT-1))
Omega <- array(NA, dim = c(3,5,nI,nT-1))

z <- matrix(NA, nI, nT)
y <- matrix(NA, nI, nT)

# simulate the population using Francis growth parameterization of the Von Bertalanffy Growth Funcation

L1 <- 160 # size (mm) at age t1 
L2 <- 403 # size (mm) at age t2
L3 <- 513 # size (mm) at age t3
rFran <- (L3 - L2) / (L2 - L1)
t1 <- 1   # age 1
t3 <- 18  # age 18
K <- (-2*log(rFran)) / (t3 - t1)
Linf <- L1 + ((L3 - L1) / (1 - (rFran)^2))

l <- s
simvals <- c(phi,k,R,p,r,L1,L2,L3)

for (i in 1:nI){
  
  s.mat[i,f[i]] <- l[i]   # assign initial size
  z[i,f[i]] <- 1 # at tagging, a fish is alive
  y[i,f[i]] <- 1 # at tagging, a fish is observed by a biologist/angler
  
  for (t in (f[i]+1):nT){
    
    s.mat[i,t] <- s.mat[i,t-1] + (K*(Linf - s.mat[i,t-1]))
    
  }
  
  for (t in f[i]:(nT-1)){  
    k.mat[i,t] <- ifelse(s.mat[i,t] >= 305, k[2], k[1]) # if big, get legal size harvest probability
    R.mat[i,t] <- ifelse(s.mat[i,t] >= 305, R[2], R[1]) # if big, get legal size catch&release prob    
    
    Psi[1:3,1:3,i,t] <- matrix(c(phi * r,     k.mat[i,t] * r,     (1 - phi - k.mat[i,t]) * r + (1 - r),
                                 0,           0,             1,
                                 0,           0,             1), nrow = 3, byrow = T)
    
    Omega[1:3,1:5,i,t] <- matrix(c(p * R.mat[i,t], p * (1 - R.mat[i,t]), (1 - p) * R.mat[i,t],  0,  (1 - p)*(1 - R.mat[i,t]),
                                   0,              0,                    0,                     1,  0,
                                   0,              0,                    0,                     0,  1), byrow = T, nrow = 3)
    
    z[i,t+1] <- which(rmultinom(1, 1, Psi[z[i,t],,i,t]) == 1)
    y[i,t+1] <- which(rmultinom(1, 1, Omega[z[i,t+1],,i,t]) == 1)
    
  }
}


# plot growth
plot(s.mat[1,] ~ seq(1,nT), type = 'l', col = 'grey50', 
     ylim = c(0,700), ylab = 'Size', xlab = 'Time')
for (i in 2:nI){
  lines(s.mat[i,f[i]:nT] ~ seq(f[i],nT), col = 'grey50')
}

# numbers by observation type and occasion
for (t in 1:nT){
  print(table(y[,t]))
}

# length data
s.dat <- matrix(NA, nI, nT)
for (i in 1:nI){
  for (t in f[i]:nT){
    if (y[i,t] != 5){
      s.dat[i,t] <- rnorm(1, s.mat[i,t], 0.1)
    }
  }
}

one <- rep(1, nI)


# Define CMR model in JAGS -------------------------------------------------------


sink("tvr_size.jags")
cat("
  model {

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # priors
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    phi ~ dbeta(1,1)
    
    R[1] ~ dbeta(1,9)
    R[2] ~ dbeta(1,9)
    
    kappa[1] ~ dbeta(1,9)
    kappa[2] ~ dbeta(1,9)
    
    p ~ dbeta(1,4)
    
    # beta priors defined by tag retention model posteriors derived outside model
    r ~ dbeta(140.2,26.8) # > rb_parms:  $alpha [1] 140.1815 $beta [1] 26.81656
    
    L1 ~ dnorm(160,10)
    L2 ~ dnorm(403,10)
    L3 ~ dnorm(513,10)


    # model growth
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 1:nI){
      s[i,f[i]] <- s.dat[i,f[i]]           # size at tagging assigned from s.dat
      for (t in (f[i] + 1):nT){
        # s[i,t] <- s[i,t-1] + g           # grow fish through time
        # s[i,t] <- s[i,t-1] + (K*(Linf - s[i,t-1]))
        s[i,t] <- s[i,t-1] + 
                  (-2 * (log((L3 - L2) / (L2 - L1))) / (t3 - t1)) *
                  ((L1 + ((L3 - L1) / (1 - ((L3 - L2) / (L2 - L1))^2))) 
                  - s[i,t-1])
        s.dat[i,t] ~ dnorm(s[i,t], 100)    # where fish size data are filled in...drawn from a normal distribution of size s at time t, where s is being estimated based on growth   
      }
      for (t in f[i]:nT){
        class[i,t] <- ifelse(s[i,t] > 305, 2, 1) # assign class based on length at time t
      }
    }  
  

    #--------------------------------------
    # Latent states (s):
    # 1) alive, with tag
    # 2) caught & kept
    # 3) unavailable
    #--------------------------------------

    for (i in 1:nI){
      for (t in f[i]:(nT-1)){
        ps[1,1,i,t] <- phi * r
        ps[1,2,i,t] <- kappa[class[i,t]] * r
        ps[1,3,i,t] <- (1 - phi - kappa[class[i,t]]) * r + (1 - r)
    
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
     }
   }

  # marginalized capture recapture likelihood
  for (i in 1:nI){      
    zeta[i,f[i],1] <- 1 
    zeta[i,f[i],2] <- 0
    zeta[i,f[i],3] <- 0

    for(t in f[i]:(nT-1)){ 
      zeta[i,(t+1),1] <- inprod(zeta[i,t,], ps[,1,i,t]) * po[1,y[i,t+1],i,t]
      zeta[i,(t+1),2] <- inprod(zeta[i,t,], ps[,2,i,t]) * po[2,y[i,t+1],i,t]
      zeta[i,(t+1),3] <- inprod(zeta[i,t,], ps[,3,i,t]) * po[3,y[i,t+1],i,t]

      }
      lik[i] <- sum(zeta[i, nT, ])
      one[i] ~ dbern(lik[i])
  }
 
}",fill = TRUE)
sink()


# Set up and run model in JAGS --------------------------------------------

# Define data
jags.data <- list(y = y, nI = nrow(y), nT = nT, f = f, one = one, s.dat = s.dat, t1 = t1, t3 = t3)

# Initial values
inits <- function(){list(phi = 0.5)}  

# Parameters monitored
parameters <- c('phi','kappa','R','p','r','s','L1','L2','L3')

# MCMC specifications
nc <- 2
nt <- 20
ni <- 5000
nb <- 1000

jags.sim <- jags(jags.data, inits, parameters, "tvr_size.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb,
                 parallel=T, n.cores = nc)
modsim <- jags.sim

print(modsim)

MCMCtrace(modsim,
          pdf = FALSE,
          params = c('phi','kappa','R','p','r','L1','L2','L3'),
          gvals = simvals,
          Rhat = TRUE,
          n.eff = TRUE,
          post_zm = TRUE)
