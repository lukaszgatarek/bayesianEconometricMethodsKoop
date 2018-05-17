#### osterior simulation in the multinomial probit, as in ex. 14.7
# there are two differences compared to the book setup. First one: we assume sigma to be known.
# Thus the Gibbs only iterates between betas and latent variable.
# On the other hand the specification of the latent process is extended in the sense that we allow for more than just
# one explanatory variable. Thus the part of Z corresponding to i-th individual is not given by

#|zi1_1 0       ...   0    |
#|0     zi1_1   ...   0    |
#|.     .    .        .    |
#|.     .     .       .    |
#|.     .      .      .    |
#|0     0       ... ziJ-1_1|

# but rather 

#|zi1_1 0       ...   0     zi1_2 0       ...   0    |
#|0     zi1_1   ...   0     0     zi1_2   ...   0    |
#|.     .    .        .     .     .    .        .    |
#|.     .     .       .     .     .     .       .    |
#|.     .      .      .     .     .      .      .    |
#|0     0       ... ziJ-1_1 0     0       ... ziJ-1_2|

# and then it trivially extends to cases more complex than 2 explanatory variables.

require(mvtnorm)
require(tmvtnorm)

# we apply a block structure for z as in albert & Chib (1993) (here z denote the exaplanatory variables in the latent process and not the depending variable as in previous latent processes cases)
# number of explanatory variables
r <- 5
# number of choices 
J <- 4
# number of individuals
n <- 100
# get the explanatory variables properly, we need r x J x n numbers for that
explanatory <- rnorm(r*J*n, 0, 1)
# prepare the container for the Z matrix; it has n x J number of rows and (r x J) columns; it is a block diagonal kind of amtrix
Z <- matrix(0, n*J, r*J)
# allocate the first explanatory variable in the matrix
for (k in 1:r){
  for (i in 1:n){
    diag(Z[( (i-1)*J + 1): (i * J),  ((k-1)*J+1):(k*J)]) <- explanatory[ ( (k-1) * (n*J) + ( (i-1)*J + 1) ) : ((k-1) * (n*J) + (i * J))]
  }
}

### simulate the error terms
# we need to get the common matrix Sigma
Sigma <- diag(J)
# we need n columns of J errors
errors <- matrix(NA, J, n)

for (i in 1:n){
  errors[,i] <- rmvnorm(1, matrix(0,J,1), Sigma)
}
# we vectorize this matrix
eps <- t(t(as.vector(errors)))

### vector beta
# we need J * r elements
trueBeta <- matrix(rnorm(J*r), J, r)

# then we need to vectorize it
beta <- t(t(as.vector(trueBeta)))

### simulate from the DGP according to eq. 14.31; the simulation is performed for all the individuals at once
# first simulate the laten process
U <- Z %*% beta +  eps
## then we convert U into y, which is the observed process
# to that end we need to get back from vectorized form of U into an individual related one (each column or one individual)
matrixU <- matrix(U, J, n)
y <- matrix(NA, n, 1)
for (i in 1:n){
  if (max(matrixU[,i])   <= 0) { y[i] <- 0 }
  if (max(0, matrixU[,i]) > 0) { y[i] <- which(matrixU[,i] == max(matrixU[,i]))}
}

# plot y against the maximum of U, for each individual
plot(apply(matrixU,2,max), y)

### Gibbs sampler (it is sampled for a model in a vectorized form)
# prior for beta
mubeta <- matrix(0, J*r, 1)
Vbeta <- diag(J*r) * 100

# containers for the draws
nSim <- 100
betaSim <- matrix(0, J * r, nSim)
USim <- matrix(0, n*J, 1)
# get some first USim to initialize the mcmc

for (i in 1:n){
    
  if (y[i] == 0){
    
    D <- matrix(0, J, J)
    diag(D) <- 1    
    mean <- as.vector(Z[((i-1)*J + 1) : (i*J), ] %*% betaSim[,1])
    upper <- rep(0,J)
    
    #    rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
    USim[ ((i-1)*J + 1) : (i*J) ] <- rtmvnorm(n=1, mean = mean, Sigma, upper = upper,  D=D, algorithm="gibbsR") 
  } else {
    
    D <- matrix(0, J, J)
    D[,y[i]] <- -1
    diag(D) <- 1
    
    mean <- as.vector(Z[((i-1)*J + 1) : (i*J), ] %*% betaSim[,1])
    upper <- rep(0,J)
    upper[y[i]] <- Inf
    USim[ ((i-1)*J + 1) : (i*J) ] <- rtmvnorm(n=1, mean = mean, Sigma, upper = upper,  D=D, algorithm="gibbsR")
  }
   
  #  plot(rtmvnorm(n=100, mean = mean, Sigma, upper = upper,  D=D, algorithm="gibbsR"))
}

# beta step
for (s in 2:nSim){
  
  print(s)
  
  # we create a huge matrix Omega, accordingly to Chib & Albert (1993)
  invOmega <-matrix(0, n * J, n * J)
  
  for (i in 1:n){
    invOmega[((i-1)*J+1):(i*J), ((i-1)*J+1):(i*J)] <- solve(Sigma) 
  }
    
  Dbeta <- solve( t(Z)%*% invOmega %*% Z + solve(Vbeta) )
  dbeta <- t(Z)%*% invOmega %*% U + solve(Vbeta) %*% mubeta
  
  betaSim[,s] <- t(rmvnorm(1, Dbeta %*% dbeta, Dbeta))
  
  # U step
  for (i in 1:n){
    
    if (y[i] == 0){
      
      D <- matrix(0, J, J)
      diag(D) <- 1    
      mean <- as.vector(Z[((i-1)*J + 1) : (i*J), ] %*% betaSim[,s])
      upper <- rep(0,J)
      
      #    rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
      USim[ ((i-1)*J + 1) : (i*J) ] <- rtmvnorm(n=1, mean = mean, Sigma, upper = upper,  D=D, algorithm="gibbsR") 
    } else {
      
      D <- matrix(0, J, J)
      D[,y[i]] <- -1
      diag(D) <- 1
      
      mean <- as.vector(Z[((i-1)*J + 1) : (i*J), ] %*% betaSim[,s])
      upper <- rep(0,J)
      upper[y[i]] <- Inf
      USim[ ((i-1)*J + 1) : (i*J) ] <- rtmvnorm(n=1, mean = mean, Sigma, upper = upper,  D=D, algorithm="gibbsR")
    }
    
    #  plot(rtmvnorm(n=100, mean = mean, Sigma, upper = upper,  D=D, algorithm="gibbsR"))
  }
}

for (m in 1:dim(beta)[1]){
  expMean <- function(j){
    return(mean(betaSim[m, 1:j]))
  }
  
  # compute the expanidng mean of alpha posterior samples
  expMeanCalc <- apply(matrix(1:dim(betaSim)[2]), 1, expMean)
  
  # plot the (exapnding) mean of the draws posterior of alpha and the true value
  plot(expMeanCalc, ylab = "expanding widow mean", xlab = "Gibbs iteration", main = paste('beta', toString(m)))
  abline(a = beta[m,1], b = 0, col = "red")
}


