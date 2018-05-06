## posterior simulation in the ordered probit model
# this is a simplified version compared to the one presented in ex 14.5
# the simplification is related to no-cutpoints step in the sampler.
# the sampler is limited to beta and z (latent variable) step, 
# what makes it very similar to the one presented in ex 14.1, 
# with the distiction that now more classes of dependent variable are concerned

rm(list = ls())

require(mvtnorm)
require(msm)

#simulate from the data generating process
T <- 1000
r <- 3 
beta <- rnorm(r, 0, 1)
# error in the latent variable equation
eps <- rnorm(T, 0, 1)
# matrix of regressors
x <- matrix(rnorm(T * r, 0, 1), T, r)
# determine the latent data
z <- x %*% t(t(beta)) + eps
## derive y from z
# first we need to select the number of cutpoints
M <- 3
a <- rep(0, M + 1) 
a0 <- -Inf
a1 <- 0
aM <- Inf
# among M cutpoints, two are always identified a prior, a1 = 0 and aM = Inf
# thus for specified M we always select only M-2 from uniform
cutpoints <- runif(M-2, 1, 5) # value of 5 is set ad hoc, can by anything reasonably bigger than 0, so one is imposed as the bottom boundary
a[1] <- -Inf 
a[2] <- 0
for (i in 3:M){
  a[i] <- cutpoints[i-2]
}
a[M + 1] <- Inf
# there is a slight notation difference between eq. 14.16 in the book and our vector. we can not have a[0] here.. thus we run up to M+1
# now we can translate our latent variable into the y's

y <- z # initialization
# if (z <= 0) {y <-1}
y[z<=0] <- 1
for (i in 2:M){
  # assign y according to the value of z
  y[(z>a[i]) & (z<a[i+1])] <- i
}

## implement the Gibbs sampler
# number of iterations in Gibbs sampler
nSim <- 1000
# prior for beta
mubeta <- matrix(0, r, 1)
Vbeta <- diag(r) * 100

# arrange variable names for the sampler
betaSim <- matrix(0, r, nSim)

# get some first z to initialize the mcmc
for (j in 1:T){
  # we sample z accordingly to the interval, where y[j] belongs to
  z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,1], sd = 1, lower = a[y[j]], upper = a[y[j] + 1])
}

for (i in 1:nSim){
 # beta step 
 Dbeta <- solve(t(x) %*% x )
 dbeta <- t(x) %*% z 
  
 # draw from a respective distribution 
 betaSim[,i] <- t(rmvnorm(1, Dbeta %*% dbeta, Dbeta))
  
 for (j in 1:T){
   
     # we sample z from the interval, which y[j] is associated to; each y is associated to another interval of sampling for z 
     z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,i], sd = 1, lower = a[y[j]], upper = a[y[j] + 1])
 }
 print(i)
}

for (i in 1:r){
  plot(density(betaSim[i,]))
  abline(v = beta[i], col = "red")
}

############### some plotting
## expanding window mean function for alphaSim

for (i in 1:r){
  expMean <- function(j){
    return(mean(betaSim[i, 1:j]))
  }
  
  # compute the expanidng mean of alpha posterior samples
  expMeanCalc <- apply(matrix(1:dim(betaSim)[2]), 1, expMean)
  
  # plot the (exapnding) mean of the draws posterior of alpha and the true value
  plot(expMeanCalc, ylab = "expanding widow mean", xlab = "Gibbs iteration", main = paste('beta', toString(i)))
  abline(a = beta[i], b = 0, col = "red")
}

