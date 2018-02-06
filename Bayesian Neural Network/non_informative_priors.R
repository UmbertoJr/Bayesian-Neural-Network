library(neuralnet)
library(matrixcalc)
library(MASS)
library(matlib)

# functions  --------------------------------------------------------------

LOG.likelihood <- function(y , Z , betas, sigma){
  n= dim(Z)[1]
  val = 0
  for(i in 1:n){
    mu <-betas[]%*%Z[i,]
    #cat(mu, '\n')
    val <- val - (y[i] - mu)^2/(2*sigma) -log(sigma)
    #cat( val,'\n')
  }
  return(val)
}

sampling_gammas_with_LL <- function(gamma.matrix, Z, y, X, betas, sigma){
  k <- dim(gamma.matrix)[1]
  p <- dim(gamma.matrix)[2]
  for(j in 1:k){
    new.G <- gamma.matrix
    new.gamma <- mvrnorm(mu = gamma.matrix[j,], Sigma = 0.05^2*diag(p))
    new.G[j,] <- new.gamma
    Z.new <- compute_Z_matrix(gamma.matrix = new.G, X.matrix = X)
    foo <- det(t(Z.new)%*%Z.new)
    #cat(c('foo'=foo, '\n'))
    #print(dim(Z.new))
    if(foo>0.001 & all(new.gamma< 1e5)){
      pr <-  LOG.likelihood(y,Z.new , betas = betas, sigma = sigma) - LOG.likelihood(y, Z , betas = betas, sigma = sigma)
      prob <- min(1,exp(pr) )
      #print(prob)
      accept <- rbinom(1,1,prob = prob)
      #print(accept)
      if( accept==1){
        Z <- Z.new
        gamma.matrix <- new.G
      }
    }
  }
  return(list('gammas' = gamma.matrix, 'Z'= Z))
}

compute_Z_matrix <- function(gamma.matrix, X.matrix){
  N = dim(X.matrix)[1]
  K = dim(gamma.matrix)[2]
  Z = matrix(nrow = N, ncol = K)
  for(i in 1:N){
    for(j in 1:K){
      Z[i,j] <- 1/ (1 + exp(gamma.matrix[1, j] + gamma.matrix[-1, j]%*% X.matrix[i,]))
    }
  }
  Z <- cbind(rep(1,N), Z)
  return(Z)
}

sampling_sigma_and_betas <- function(Z,y){
  # generating sigma squared
  n <- dim(Z)[1]
  k <- dim(Z)[2]
  alpha <- (n-k-1)/2
  beta <- 1/2 * t(y)%*%(diag(n)- Z%*%inv(t(Z)%*%Z)%*%t(Z))%*%y
  tau.2 <- rgamma(1, alpha, rate = beta)
  
  # generating Betas
  foo <- inv(t(Z)%*%Z)
  Mu <- foo%*%t(Z)%*%y
  Sigma <- foo / tau.2
  Betas <- mvrnorm(mu = Mu, Sigma = Sigma)
  return(list('sigma'= 1/sqrt(tau.2) , 'Betas'=Betas))
}


# simulation --------------------------------------------------------------

# Dataset creation
# N=100     # number of observations
# P=5       # number of covariates
# K = 20    # number of neurons in the hidden layer
# 
# X <- matrix(runif(N*P), c(N,P))
# Y <- (X[,1]^2 + log(X[,2]) + 3^X[,3]+ X[,4]*X[,5]) + c(rnorm(100))
# dat <- data.frame(cbind(X,Y))
#save(file = 'data_simmulate.RData', dat, X,Y)

load(file = 'data_simmulate.RData')

#classical NN to initialize
K=18
f <- as.formula(paste("Y~" , paste(names(dat)[!names(dat)=='Y'], collapse = '+')))
net_for_conj_priors <- neuralnet(f, data = dat, hidden = K , linear.output = T,
                   lifesign = 'full', threshold = 0.1)

plot(net_for_conj_priors$net.result[[1]], Y)
sum((net_for_conj_priors$net.result[[1]]- Y)^2)/42

# MCMC
p = dim(X)[2] + 1
X.matrix<- as.matrix(X)
y.vector <- as.array(Y)
gamma.matrix <- net_for_conj_priors$weights[[1]][[1]]
Z <- compute_Z_matrix(gamma.matrix, X.matrix)
det(t(Z)%*%Z)
conjugate_sampling <- sampling_sigma_and_betas(Z,y.vector) 
Betas <- conjugate_sampling$Betas
sigma <- conjugate_sampling$sigma

LOG.likelihood(y=Y,Z=Z, betas = Betas, sigma=sigma)

MH_sampling <- sampling_gammas_with_LL(gamma.matrix = gamma.matrix,
                                       Z = Z,y = y.vector,X = X.matrix, betas = Betas, sigma = sigma)



simualtion_BNN_flat_prior <- function(X, y,K, Betas, Gammas, sigma, 
                                      it = 1e3, burn.in=1e3, tm = 3, plot.at=100){
  P= dim(X)[2]
  chain <- data.frame(matrix(NA, nrow = it,ncol = (K*(P+1)+ K+1 +1)))
  Z <- compute_Z_matrix(gamma.matrix= Gammas, X.matrix=X)
  t = 1
  h = 1
  pl = 1
  burn = 0
  sim <- it*tm + burn.in
  print(sim)
  for( i in 1:sim){
    conjugate_sampling <- sampling_sigma_and_betas(Z,y.vector) 
    Betas <- conjugate_sampling$Betas
    sigma <- conjugate_sampling$sigma
    MH_sampling <- sampling_gammas_with_LL(gamma.matrix = Gammas,
                                           Z = Z,y = y.vector,X = X.matrix,
                                           betas = Betas, sigma = sigma)
    
    Z <- MH_sampling$Z
    Gammas <- MH_sampling$gammas
    if(i > burn.in){
      burn <- 1
      t = t + 1
    }
    if(((t==tm) && (burn == 1))){
      chain[h,1:(K+1)] <- conjugate_sampling$Betas
      chain[h, (K*(P+1)+ K+1 +1)] <- conjugate_sampling$sigma
      chain[h,(K+2):(K+1+(K*(P+1)))] <- as.vector(MH_sampling$gammas)
      t = 0
      if(pl == plot.at){
        print(h)
        pl = 0
      }
      h = h +1
      pl = pl + 1
    }
  }
  return(chain)
}

names <- rep(NA, (K*(P+1)+ K+1 +1))
for(l in 0:(K+1)){
  print(l)
  names[l+1] <- paste('Beta',l)
}
l=K+2
for(i in 1:K){
  for(p in 1:(P+1)){
    names[l] <- paste('Gamma',paste(i,p,sep = '.'))
    l = l + 1
  }
}
names[(K*(P+1)+ K+1 +1)] <- 'sigma'




c1 <- simualtion_BNN_flat_prior(X,Y,K=K, Betas, Gammas = gamma.matrix,
                                sigma = sigma ,it = 10000 ,
                                burn.in = 1000 ,tm = 5, plot.at=100)

colnames(c1) <- names
#save(c1, file = 'non_informative_chain.RData')
load('non_informative_chain.RData')




# functions to predict ----------------------------------------------------

predicting_with_join <- function(X,y, chain){
  lung <- dim(chain)[1]
  for(i in 1:lung){
    betas <- chain[i,seq(from = 1,to = 19)]
    Gammas <- chain[i,seq(from = 20,to = 19)]
    val <- LOG.likelihood(y,Z,betas,sigma)
  }
}


