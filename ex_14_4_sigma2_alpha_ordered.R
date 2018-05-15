## implementation of exercise 14.4 from Koop et al. Bayesian Econometric Methods (with general sigma2_alpha)
# posterior simulation in the panel probit model, with general sigma_alpha;
# the difference between this implementation and the previous one (ex_14_4) results from additional parameter, 
# sigma2_alpha, which enters the equation describing alpha i.e. alpha_i ~ N(alpha, sigma2_alpha), 
# but that has been already implemented in ex_14_4_sigma2_alpha

# Here, we add one more complication in the sense that the model is not a probit, but rather an ordered probit model, 
# which makes it difficult to be estimated. This difficulty follows from dense information on latent variable z. 
# as the domain is cut with cutpoints, in each interation we get extremely precise coverage of each interval between two cutpoints.
# that makes it almost impossible to move a cutpoint in an subsequent iteration, because the sampler uses formulas

    ###   alpha_m <- max(cutpoints[m - 1, s - 1], max(z[y == m - 1]))
    ###   beta_m <- min(cutpoints[m + 1, s - 1], min(z[y == m])) 
    ###   cutpoints[m, s] <- runif(1, alpha_m, beta_m)

# where again z's in the current iteration were based on cutpoints[, s-1], and thus the current draw of cutpoints is close to the previous one and so on.
# to milden this problem, we have decided not to use z entirely, but rather to use a sample from it. This way we allow for some randomness in cutpoint level, 
# as some z's close to the cutpoints are just not taken into account, and in consequence, the calculated min or max of z's do not reflect its closeness to the cutpoint,
# allowing for some improved draw of a cutpoint in a given iteration

rm(list = ls())

require(mvtnorm)
require(msm)
require("PerformanceAnalytics")
require('invgamma')

# simulate from the data generating process
T <- 100
# number of explanatory variables in latent variable equation (common for every i); those are variables corresponding only to beta (alpha is not included in r)
r <- 2 
## select the number of individuals in the panel
# the accuracy of estimation grows with n; for low n, the accuracy is usually very poor in terms of estimating alpha parameter
n <- 50 
# select beta parameter that is common across individuals
beta <- rnorm(r, 0, 1)
## select alphas which is associated with i-th individual 
# to that end we need some value of alpha, which drives the distribution of alpha_i's (all alpha_i's are centered around alpha)
alpha <- rnorm(1, 0, 0.5)
# in the general model sigma_alpha_2 is also a model parameter, we draw its value (for DGP) from inversegamma distribution
sigma2_alpha <- rinvgamma(1,150,100)
# get some values of alpha_i for each individual
alpha_i <- rnorm(n, alpha, sqrt(sigma2_alpha)) # for the time being sigma2_alpha is fixed at 1
# to work with the stacked, panel data structure, we need to repeat the values of each alpha_i element T times
alpha_is <- matrix(t(matrix(alpha_i, length(alpha_i), T)))
# error in the latent variable equation (stacked i-th by (i-1)-th, forming a Txn length vector)
eps <- as.matrix(rnorm(n * T, 0, 1))
# matrix of regressors (again, we stack them indivdidual by individual)
X <- matrix(rnorm(n * T * r, 0, 0.1), n * T, r) # we simulate n x T x r values and then we structure them into nT x r matrix
# having all the components, we can simulate from this data generating process according to 14.10
z <- alpha_is + X %*% beta + eps
## derive y from z
# first we need to select the number of cutpoints
M <- 7
a <- rep(0, M + 1) 
a0 <- -Inf
a1 <- 0
aM <- Inf
# among M cutpoints, two are always identified a prior, a1 = 0 and aM = Inf
# thus for specified M we always select only M-2 from uniform
trueCutpoints <- runif(M-2, 0.1, 2) # value of 2 is set ad hoc, can by anything reasonably bigger than 0
# reorder in an increasing way
trueCutpoints <- trueCutpoints[order(trueCutpoints)]

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
mubeta <- matrix(0, r, 1) # there is r values in beta vector, each corresponding to different variable x
Vbeta <- diag(r) * 100 # this is a variance beta matrix; following example in ex. 14.2 we assume diag(r) x 100
# prior for alpha
mualpha <- 0
Valpha <- 1
# prior for sigma2
a_prior <- 3
b <- .5

# arrange variables for the sampler
betaSim <- matrix(0, r, nSim) # all individuals have the same parameter vector beta
alphaiSim <- matrix(0, n, nSim) # every individual has different value
alphaSim <- matrix(0, 1, nSim) # there is a common value which drives all alpha_i's
sigmasAlpha2Sim <- matrix(1, 1, nSim) # apart from common mean, alpha, there is parameter of variance, sigmaAlpha2, which drives the dispersion of the individual alphas

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

# # get some first z to initialize the mcmc (we need some initial value of cutpoints to initialize the latent variables)
# for (j in 1:T){
#   # we sample z accordingly to the truncated interval, y[j] is associated with
#   z[j] <- rtnorm(1, mean = x[j,] %*% betaSim[,1], sd = 1, lower = cutpoints[y[j], 1], upper = cutpoints[y[j] + 1, 1])
# }

## we need to start the Gibbs somewhere, let us  get some first z's to initialize the mcmc
# formula 14.11
for(i in 1:n){  # for each individual
  for (j in 1:T){ # for each time series dimension observation
    observation <- ((i-1) * T) + j # values are stack one on another [T' T' ... T']'
    # if (y[observation] == 0){ 
    #   # sample from the truncated normal [-Inf, 0]
    #   z[observation] <- rtnorm(1, mean = alphaiSim[i,1] + X[observation,] %*% betaSim[,1], sd = 1, lower = -Inf, upper = 0)
    #       } else if (y[observation] == 1){
    #   # sample from the truncated normal [0, Inf]     
    #   z[observation] <- rtnorm(1, mean = alphaiSim[i,1] + X[observation,] %*% betaSim[,1], sd = 1, lower = 0,    upper = Inf)
    #       }
    # 
    z[observation] <- rtnorm(1, mean = alphaiSim[i,1] + X[observation,] %*% betaSim[,1], sd = 1, lower = cutpoints[y[observation], 1], upper = cutpoints[y[observation] + 1, 1])
  }
}

## after getting some initial vector of z's, we can run the mcmc sampler
for (s in 2:nSim){
    
  print(s)
  
   ## alpha_i's step according to the formula 14.11 and 14.12
   for (i in 1:n){
      
      D_alpha_i <- (T + (sigmasAlpha2Sim[s-1])^(-1) )^(-1)
      d_alpha_i <- sum(z[ ((i-1) * T + 1) : (i * T) ] - X[((i-1) * T + 1) : (i * T), ] %*% betaSim[, s-1]) + (sigmasAlpha2Sim[s-1])^(-1) * alphaSim[s-1]
      
      alphaiSim[i, s] <- rnorm(1, D_alpha_i %*% d_alpha_i, sqrt(D_alpha_i) )      
   }    
    
   ## beta step according to the formula 14.13 
   alpha_is <- matrix(t(matrix(alphaiSim[, s], length(alphaiSim[, s]), T)))
    
   Dbeta <- solve(t(X) %*% X + solve(Vbeta) )
   dbeta <- t(X) %*% (z - alpha_is) + solve(Vbeta) %*% mubeta
    
   # draw from a respective distribution 
   betaSim[,s] <- t(rmvnorm(1, Dbeta %*% dbeta, Dbeta) )
    
   ## alpha step according to the fomrula 14.14
   D_alpha <- (n * (sigmasAlpha2Sim[s-1])^(-1) + (Valpha)^(-1) )^(-1)
   d_alpha <- sum(alphaiSim[, s]) * (sigmasAlpha2Sim[s-1])^(-1) + Valpha^(-1) * mualpha 
   
   alphaSim[s] <- rnorm(1, D_alpha * d_alpha, sqrt(D_alpha))

   ## sigma_alpha step according to 14.15
   sigmasAlpha2Sim[s] <- rinvgamma(1, n/2 + a_prior, 1/b + .5 * t(alphaiSim[, s] - alphaSim[s]) %*% (alphaiSim[, s] - alphaSim[s]) )
   
   ## z step accroding to the formula 14.11
   
   for(i in 1:n){  # for each individual
     for (j in 1:T){
       observation <- ((i-1) * T) + j
       
       # if (y[observation] == 0){ 
       #   #    rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
       #   z[observation] <- rtnorm(1, mean = alphaiSim[i,s] + X[observation,] %*% betaSim[,s], sd = 1, lower = -Inf, upper = 0)
       # } else if (y[observation] == 1){
       #   z[observation] <- rtnorm(1, mean = alphaiSim[i,s] + X[observation,] %*% betaSim[,s], sd = 1, lower = 0,    upper = Inf)
       # }
       
       z[observation] <- rtnorm(1, mean = alphaiSim[i,s] + X[observation,] %*% betaSim[,s], sd = 1, lower = cutpoints[y[observation], s-1], upper = cutpoints[y[observation] + 1, s-1])
       
     }
   }
   
   # the convergence is to slow due precision z covers the intervals between the cutpoints.
   # As it covers exactly the intervals between the cutpoints, it is unlikely that the new 
   # drawn cutpoint gets changed significantly from iteration to the iteration.
   # Therefore instead of sampling cutpoints based on full z, we just take a sample from z's for inference on the new cutpoint
   
   indicators <- sample(1:length(z), length(z)/10)
   zSampled <- z[indicators]
   ySampled <- y[indicators]
   
   for (m in 3:M){ # in fact this m here shall go from 2:M, but we can not have the index of 0 i.e. cutpoints[0,]
     # thus, we have m - 1 and m instead of m and m + 1
 
     # alpha_m <- max(cutpoints[m - 1, s - 1], max(z[y == m - 1]))
     # beta_m <- min(cutpoints[m + 1, s - 1], min(z[y == m])) 
         
     alpha_m <- max(cutpoints[m - 1, s - 1], max(zSampled[ySampled == m - 1]))
     beta_m <- min(cutpoints[m + 1, s - 1], min(zSampled[ySampled == m])) 
     cutpoints[m, s] <- runif(1, alpha_m, beta_m)
   }
   
}

## plot the densities of betas' samples vs the true value
for (i in 1:r){
  plot(density(betaSim[i,]), main = paste('beta_', toString(i)))
  abline(v = beta[i], col = "red")
}

## plot the densities of individual alphas samples vs the true value
for (i in 1:n){
  plot(density(alphaiSim[i,]), main = paste('alpha_', toString(i)))
  abline(v = alpha_i[i], col = "red")
}


############### some plotting
## expanding window mean function for alphaSim
expMean <- function(i){
  return(mean(alphaSim[1,1:i]))
}

# compute the expanidng mean of alpha posterior samples
expMeanCalc <- apply(matrix(1:dim(alphaSim)[2]), 1, expMean)

# plot the (exapnding) mean of the draws posterior of alpha and the true value
plot(expMeanCalc)
abline(a = alpha, b = 0, col = "red")


## expanding window mean function for simgasAlpha2Sim
expMean <- function(i){
  return(mean(sigmasAlpha2Sim[1,1:i]))
}

# compute the expanidng mean of alpha posterior samples
expMeanCalc <- apply(matrix(1:dim(sigmasAlpha2Sim)[2]), 1, expMean)

# plot the (exapnding) mean of the draws posterior of sigmasAlpha and the true value
plot(expMeanCalc)
abline(a = sigma2_alpha, b = 0, col = "red")

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



