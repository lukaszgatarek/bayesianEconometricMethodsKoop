## posterior simulation in the probit model
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
x <- matrix(rnorm(T * r, 0, 10), T, r)
# determine the latent data
z <- x %*% t(t(beta)) + eps
# derive y from z
y <- (z > 0)
y[y==TRUE] <- 1
y[y==FALSE] <- 0

## implement the Gibbs sampler
# number of iterations in Gibbs sampler
nSim <- 100
# prior for beta
mubeta <- matrix(0, r, 1)
Vbeta <- diag(r) * 100

# arrange variable names for the sampler
betaSim <- matrix(0, r, nSim)

# get some first z to initialize the mcmc
for (j in 1:T){
  
  if (y[j] == 0){ 
#    rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
    z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,i], sd = 1, lower = -Inf, upper = 0)
  } else if (y[j] == 1){
    z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,i], sd = 1, lower = 0,    upper = Inf)
  }
}

for (i in 1:nSim){
 # beta step 
 Dbeta <- solve(t(x) %*% x + solve(Vbeta) )
 dbeta <- t(x) %*% z + solve(Vbeta) %*% mubeta
  
 # draw from a respective distribution 
 betaSim[,i] <- t(rmvnorm(1, Dbeta %*% dbeta, Dbeta))
  
 for (j in 1:T){
   
    if (y[j] == 0){ 
      z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,i], sd = 1, lower = -Inf, upper = 0)
    } else if (y[j] == 1){
      z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,i], sd = 1, lower = 0,    upper = Inf)
    }
 }
 print(i)
}

for (i in 1:r){
  plot(density(betaSim[i,]))
  abline(v = beta[i], col = "red")
}


