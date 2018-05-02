# implementation of exercise 14.4 from Koop et al. Bayesian Econometric Methods
## posterior simulation in the panel probit model
rm(list = ls())

require(mvtnorm)
require(msm)
require("PerformanceAnalytics")

# simulate from the data generating process
T <- 100
# number of explanatory variables in latent variable equation (common for every i); those are variables corresponding only to beta (alpha is not included in r)
r <- 3 
## select the number of individuals in the panel
# the accuracy of estimation grows with n; for low n, the accuracy is usually very poor in terms of estimating alpha parameter
n <- 200 
# select beta parameter that is common across individuals
beta <- rnorm(r, 0, 1)
## select alphas which is associated with i-th individual 
# to that end we need some value of alpha, which drives the distribution of alpha_i's (all alpha_i's are centered around alpha)
alpha <- rnorm(1, 0, 0.5)
#alpha <- 0.4
# in the general model sigma_alpha_2 is also a model parameter, here we just assume that it is equal to 1
sigma2_alpha <- 1
alpha_i <- rnorm(n, alpha, sigma2_alpha) # for the time being sigma2_alpha is fixed at 1
# to work with the stacked structure, we need to repeat the values of each alpha_i element T times
alpha_is <- matrix(t(matrix(alpha_i, length(alpha_i), T)))
# error in the latent variable equation (stacked i-th by (i-1)-th, forming a Txn length vector)
eps <- as.matrix(rnorm(n * T, 0, 1))
# matrix of regressors (again, we stack them indivdidual by individual)
X <- matrix(rnorm(n * T * r, 0, 0.1), n * T, r) # we simulate n x T x r values and then we structure them into nT x r matrix
# having all the components, we can simulate from this model
z <- alpha_is + X %*% beta + eps
# derive y from z
y <- (z > 0)
y[y==TRUE]  <- 1
y[y==FALSE] <- 0
## implement the Gibbs sampler
# number of iterations in Gibbs sampler
nSim <- 100
# prior for beta
mubeta <- matrix(0, r, 1) # there is r values of beta, each corresponding to a distinct variable x
Vbeta <- diag(r) * 100 # this is a variance beta matrix; following example in ex. 14.2 we assume diag(r) x 100
# prior for alpha
mualpha <- 0
Valpha <- 1

# arrange variables for the sampler
betaSim <- matrix(0, r, nSim) # all individuals have the same parameter vector beta
alphaiSim <- matrix(0, n, nSim) # every individual has different value
alphaSim <- matrix(0, 1, nSim) # there is a common value which drives all alpha_i's

## get some first z's to initialize the mcmc
for(i in 1:n){  # for each individual
  for (j in 1:T){
    observation <- ((i-1) * T) + j # values are stack one on another [T' T' ... T']'
    if (y[observation] == 0){ 
  #    rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
      z[observation] <- rtnorm(1, mean = alphaiSim[i,1] + X[observation,] %*% betaSim[,1], sd = 1, lower = -Inf, upper = 0)
    } else if (y[observation] == 1){
      z[observation] <- rtnorm(1, mean = alphaiSim[i,1] + X[observation,] %*% betaSim[,1], sd = 1, lower = 0,    upper = Inf)
    }
  }
}

# after getting some initial vector of z's, we can run the mcmc sampler
for (s in 2:nSim){
  
 # alpha_i's step
 for(i in 1:n){
    D_alpha_i <- (T + 1)^(-1)
    d_alpha_i <- sum(z[ ((i-1) * T + 1) : (i * T) ] - X[((i-1) * T + 1) : (i * T), ] %*% betaSim[, s-1]) + alphaSim[s-1]
    
    alphaiSim[i, s] <- rnorm(1, D_alpha_i %*% d_alpha_i, sqrt(D_alpha_i) )      
 }    
  
 # beta step 
 alpha_is <- matrix(t(matrix(alphaiSim[, s], length(alphaiSim[, s]), T)))
  
 Dbeta <- solve(t(X) %*% X + solve(Vbeta) )
 dbeta <- t(X) %*% (z - alpha_is) + solve(Vbeta) %*% mubeta
  
 # draw from a respective distribution 
 betaSim[,s] <- t(rmvnorm(1, Dbeta %*% dbeta, Dbeta) )
  
 # alpha step
 D_alpha <- (n + (Valpha)^(-1) )^(-1)
 d_alpha <- sum(alphaiSim[, s]) + Valpha^(-1) * mualpha 
 
 alphaSim[s] <- rnorm(1, D_alpha * d_alpha, sqrt(D_alpha))
 # for (j in 1:T){
 #   
 #    if (y[j] == 0){ 
 #      z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,s], sd = 1, lower = -Inf, upper = 0)
 #    } else if (y[j] == 1){
 #      z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,s], sd = 1, lower = 0,    upper = Inf)
 #    }
 # }
 
 for(i in 1:n){  # for each individual
   for (j in 1:T){
     observation <- ((i-1) * T) + j
     if (y[observation] == 0){ 
       #    rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
       z[observation] <- rtnorm(1, mean = alphaiSim[i,s] + X[observation,] %*% betaSim[,s], sd = 1, lower = -Inf, upper = 0)
     } else if (y[observation] == 1){
       z[observation] <- rtnorm(1, mean = alphaiSim[i,s] + X[observation,] %*% betaSim[,s], sd = 1, lower = 0,    upper = Inf)
     }
   }
 }
 
 print(s)
}

for (i in 1:r){
  plot(density(betaSim[i,]), main = paste('beta_', toString(i)))
  abline(v = beta[i], col = "red")
}

for (i in 1:n){
  plot(density(alphaiSim[i,]), main = paste('alpha_', toString(i)))
  abline(v = alpha_i[i], col = "red")
}

apply(alphaiSim[,500:nSim], 1, mean)

mean(alphaSim)


expMean <- function(i){
  return(mean(alphaSim[1,1:i]))
}
expMeanCalc <- apply(matrix(1:dim(alphaSim)[2]), 1, expMean)
plot(expMeanCalc)
abline(a = alpha, b = 0, col = "red")


