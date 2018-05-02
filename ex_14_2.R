## posterior simulation in the probit model ex.14.2
rm(list = ls())

require(mvtnorm)
require(msm)
require('invgamma')

#simulate from the data generating process
T <- 1000
# select the number of explanatory variables
r <- 3 
beta <- rnorm(r, 0, 1)
sigma2 <- 0.3 
# error in the latent variable equation
eps <- rnorm(T, 0, sqrt(sigma2))
# matrix of regressors
x <- matrix(rnorm(T * r, 0, 1), T, r)
# determine the latent data
z <- x %*% t(t(beta)) + eps
# derive y from z
y <- (z > 0)
y[y==TRUE] <- 1
y[y==FALSE] <- 0

## implement the Gibbs sampler
# number of iterations in Gibbs sampler
nSim <- 1000
# prior for beta
mubeta <- matrix(0, r, 1)
Vbeta <- diag(r) * 100
# prior for sigma2
a <- 3
b <- .5

# arrange variable names for the sampler
betaSim <- matrix(0, r, nSim)
sigmas2Sim <- matrix(1, nSim, 1)

# get some first z to initialize the mcmc
for (j in 1:T){
 # despite the sigma presence in the eps in the z equation, we initialize with the standard deviation of 1 
  if (y[j] == 0){ 
#    rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
    z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,1], sd = 1, lower = -Inf, upper = 0)
  } else if (y[j] == 1){
    z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,1], sd = 1, lower = 0,    upper = Inf)
  }
}

for (i in 2:nSim){
 # beta step 
 Dbeta <- solve(t(x) %*% x / sigmas2Sim[i-1] + solve(Vbeta) )
 dbeta <- t(x) %*% z / sigmas2Sim[i-1] + solve(Vbeta) %*% mubeta
  
 # draw from the respective Normal distribution 
 betaSim[,i] <- t(rmvnorm(1, Dbeta %*% dbeta, Dbeta))
  
 for (j in 1:T){
   
    if (y[j] == 0){ 
      z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,i], sd = sqrt(sigmas2Sim[i-1]), lower = -Inf, upper = 0)
    } else if (y[j] == 1){
      z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,i], sd = sqrt(sigmas2Sim[i-1]), lower = 0,    upper = Inf)
    }
 }
 
 sigmas2Sim[i] <- rinvgamma(1, T/2 + a, 1/b + .5 * t(z - x %*% betaSim[,i]) %*% (z - x %*% betaSim[,i]) )
 
 print(i)
}

## plot the draws of each beta and 
# apparently the samples frm posterior tend to be multimodally distributed; that should be researched further. 
# still the mean of the sample is rather well estimating the true value of the parameters driving the data generating process
for (i in 1:r){
  plot(density(betaSim[i,]))
  abline(v = beta[i], col = "red")
}

# get the mean of posterior draws for each beta
apply(betaSim, 1, mean)


## identification problem
# the further research shall concern the issue of identification. normally we shall observe a slight identification problem
# while estimating beta and sigma at the same time. this is a theoretical issue. in simulations with selected parameters almost unobserved
# still, the chain can exhibit some patterns of instability and it is definitely convergent if we work with a chain of beta/sigma draws, instead
# TO DO: research the issue of accuracy of beta chain vs beta/sigma chain

#mean( betaSim[1,2:nSim] / (sqrt(sigmas2Sim[2:nSim])) )
#mean( betaSim[2,2:nSim] / (sqrt(sigmas2Sim[2:nSim])) )
#mean( betaSim[3,2:nSim] / (sqrt(sigmas2Sim[2:nSim])) )