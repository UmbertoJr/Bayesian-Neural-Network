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
    val <- val - (y[i] - mu)^2/(2*sigma) -1/2 * log(sigma^2)
    #cat( val,'\n')
  }
  return(val)
}

sampling_gammas_with_LL <- function(gamma.matrix, Z, y, X, betas, sigma, coeff =0.05){
  k <- dim(gamma.matrix)[1]
  p <- dim(gamma.matrix)[2]
  for(j in 1:k){
    new.G <- gamma.matrix
    new.gamma <- mvrnorm(mu = gamma.matrix[j,], Sigma = coeff^2*diag(p))
    new.G[j,] <- new.gamma
    Z.new <- compute_Z_matrix(gamma.matrix = new.G, X.matrix = X)
    foo <- det(t(Z.new)%*%Z.new)
    #cat(c('foo'=foo, '\n'))
    #print(dim(Z.new))
    if(foo>0.001 & all(new.gamma< 1e5)){
      pr <- LOG.likelihood(y,Z.new , betas = betas, sigma = sigma) - 
        LOG.likelihood(y, Z , betas = betas, sigma = sigma)
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

sampling_sigma_and_betas_with_priors <- function(Z,y, mu_beta, sigma_beta, a.0, b.0){  #modifica con priors
  # setting hyper-parameters
  n <- dim(Z)[1]
  k <- dim(Z)[2]
  alpha <- a.0 + n/2
  foo <- (t(Z)%*%Z+ sigma_beta)
  foo.i <- inv(foo)
  Mu <- foo.i%*%(t(Z)%*%y + sigma_beta%*%mu_beta)
  
  beta <- b.0 + 1/2 *( t(y)%*%y + t(mu_beta)%*%sigma_beta%*%mu_beta - t(Mu)%*%foo%*%Mu)
  
  # generating sigma squared
  tau.2 <- rgamma(1, alpha, rate = beta)
  
  Sigma <- foo.i / tau.2
  # generating Betas
  Betas <- mvrnorm(mu = Mu, Sigma = Sigma)
  return(list('sigma'= 1/sqrt(tau.2) , 'Betas'=Betas))
}


# p = dim(X)[2] + 1
# X.matrix<- as.matrix(X)
# y.vector <- as.array(Y)
# gamma.matrix <- net_for_conj_priors$weights[[1]][[1]]
# #gamma.matrix[1,] <- net_for_conj_priors$weights[[1]][[1]][6,] # putting the bias in the first row
# #gamma.matrix[5,] <- net_for_conj_priors$weights[[1]][[1]][1,]
# mu_beta <- net_for_conj_priors$weights[[1]][[2]]
# sigma_beta <- 0.05^2*diag(K+1)
# 
# Z <- compute_Z_matrix(gamma.matrix, X.matrix)
# det(t(Z)%*%Z)
# conjugate_sampling <- sampling_sigma_and_betas_with_priors(Z,y.vector,mu_beta,sigma_beta,a.0 = 41, b.0 = 0.12) 
# Betas <- conjugate_sampling$Betas
# sigma <- conjugate_sampling$sigma
# 
# LOG.likelihood(y=Y,Z=Z, betas = Betas, sigma=sigma)
# 
# MH_sampling <- sampling_gammas_with_LL(gamma.matrix = gamma.matrix,
#                                        Z = Z,y = y.vector,X = X.matrix, betas = Betas, sigma = sigma)

simualtion_BNN_with_pr <- function(X, y,K, Betas, Gammas, sigma, 
                                      it = 1e3, burn.in=1e3, tm = 3, plot.at=100){
  P= dim(X)[2]
  no= dim(X)[1]
  mu_beta <- Betas
  sigma_beta <- 0.05^2*diag(K+1)
  chain <- data.frame(matrix(NA, nrow = it,ncol = (K*(P+1)+ K+1 +1)))
  Y.rep <- data.frame(matrix(NA, nrow = it,ncol = no))
  dev <- rep(NA, it)
  Z <- compute_Z_matrix(gamma.matrix= Gammas, X.matrix=X)
  t = 1
  h = 1
  pl = 1
  burn = 0
  sim <- it*tm + burn.in
  print(sim)
  for( i in 1:sim){
    conjugate_sampling <- sampling_sigma_and_betas_with_priors(Z,y.vector,
                                                               mu_beta,sigma_beta,
                                                               a.0 = 41, b.0 = 0.12) 
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
      p.log <- rep(NA, no)
      for(el in 1:no){
        mu <- Z[el,]%*%conjugate_sampling$Betas
        Y.rep[h,el] <- rnorm(n = 1,mean = mu, sd = conjugate_sampling$sigma)
        p.log[el] <- - (y[el] - mu)^2/(2*conjugate_sampling$sigma) -1/2*log(conjugate_sampling$sigma^2)
      }
      dev[h] <- -2*sum(p.log)
      t = 0
      if(pl == plot.at){
        print(h)
        pl = 0
      }
      h = h +1
      pl = pl + 1
    }
  }
  return(list("chain"=chain, "Y.rep"=Y.rep, "deviance"=dev))
}






# Simulation --------------------------------------------------------------

# load data
load(file = 'data_simmulate.RData')
K=18
f <- as.formula(paste("Y~" , paste(names(dat)[!names(dat)=='Y'], collapse = '+')))
net_for_conj_priors <- neuralnet(f, data = dat, hidden = K , linear.output = T,
                   lifesign = 'full', threshold = 0.1)


c3 <- simualtion_BNN_with_pr(X,Y,K=K, Betas, Gammas = gamma.matrix,
                          sigma = sigma ,it = 10000 ,
                          burn.in = 1000 ,tm = 5, plot.at=100)

colnames(c3$chain) <- names
#save(c3, file = 'chain_conjugate_BNN_2.RData')

#load('chain_conjugate_BNN.RData')
#load('chain_conjugate_BNN_2.RData')

plot(density(c2$`Beta 0`))
net_for_conj_priors$weights[[1]][[2]]
