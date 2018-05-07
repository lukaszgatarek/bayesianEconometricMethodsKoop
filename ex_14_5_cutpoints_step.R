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
M <- 4
a <- rep(0, M + 1) 
a0 <- -Inf
a1 <- 0
aM <- Inf
# among M cutpoints, two are always identified a prior, a1 = 0 and aM = Inf
# thus for specified M we always select only M-2 from uniform
trueCutpoints <- runif(M-2, 0.1, 2) # value of 2 is set ad hoc, can by anything reasonably bigger than 0

a[1] <- -Inf 
a[2] <- 0
for (m in 3:M){
  a[m] <- trueCutpoints[m-2]
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

## arrange variable names for the sampler
# for beta parameters random draws 
betaSim <- matrix(0, r, nSim)
# for cutpoints random draws
cutpoints <- matrix(0, M + 1, nSim)
# we need to insert a0, a1 and a(M+1) as they are invariant
cutpoints[1,] <- -Inf
cutpoints[2,] <- 0
cutpoints[M+1,] <- Inf

for (m in 3:M){
  # we just get some random value from uniform (1,2) to start the chain somewhere
  cutpoints[m, 1] <- runif(1,1,2)
}

# if the chain starting values are not ordered properly, then we reorder them
cutpoints[, 1] <- cutpoints[order(cutpoints[, 1]), 1] 

# get some first z to initialize the mcmc (we need some initial value of cutpoints to initialize the latent variables)
for (j in 1:T){
  # we sample z accordingly to the truncated interval, y[j] is associated with
  z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,1], sd = 1, lower = cutpoints[y[j], 1], upper = cutpoints[y[j] + 1, 1])
}

for (i in 2:nSim){
  # beta step 
  Dbeta <- solve(t(x) %*% x )
  dbeta <- t(x) %*% z 
  
  # draw from a respective distribution 
  betaSim[,i] <- t(rmvnorm(1, Dbeta %*% dbeta, Dbeta))
  
  # z step
  for (j in 1:T){
    
    # we sample z from the interval, which y[j] is associated to; each y is associated to another interval of sampling for z 
    z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,1], sd = 1, lower = cutpoints[y[j], i-1], upper = cutpoints[y[j] + 1, i-1])
  }
  
  # cutpoints step
  
  
  for (m in 3:M){ # in fact this m here shall go from 2:M, but we can not have the index of 0 i.e. cutpoints[0,]
                  # thus, we have m - 1 and m instead of m and m + 1
 #   alpha_m <- max(cutpoints[m - 1, i], max(z[y == m - 1]))
 #   beta_m <- min(cutpoints[m + 1, i], min(z[y == m]))
    
    alpha_m <- max(cutpoints[m - 1, i], max(z[y == m - 1]))
    beta_m <- min(cutpoints[m + 1, i-1], min(z[y == m])) # cutpoints[m + 1, i ] has not been available yet at computing beta_m, so we take the one from i-1 th draw
    cutpoints[m, i] <- runif(1, alpha_m, beta_m)
  }
  
  if (i %% 100 == 0 ) {print(i)}
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


for (m in 3:M){
  expMean <- function(j){
    return(mean(cutpoints[m, 1:j]))
  }
  
  # compute the expanidng mean of alpha posterior samples
  expMeanCalc <- apply(matrix(1:dim(cutpoints)[2]), 1, expMean)
  
  # plot the (exapnding) mean of the draws posterior of alpha and the true value
  plot(expMeanCalc, ylab = "expanding widow mean", xlab = "Gibbs iteration", main = paste('cutpoint', toString(m-1)))
  abline(a = a[m], b = 0, col = "red")
}


