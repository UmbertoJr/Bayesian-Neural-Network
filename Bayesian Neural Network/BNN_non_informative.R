
library(matlib)
library(MASS)
library(matrixcalc)
library(neuralnet)

# Building functions ------------------------------------------------------------------------

compute_Z_matrix <- function(gamma.matrix, X.matrix){
  N = dim(X.matrix)[1]
  K = dim(gamma.matrix)[1]
  Z = matrix(nrow = N, ncol = K)
  for(i in 1:N){
    for(j in 1:K){
      Z[i,j] <- 1/ (1 + exp(gamma.matrix[j,1] + gamma.matrix[j,-1]%*% X.matrix[i,]))
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
  sigma.2 <- rgamma(1, alpha, rate = beta)
  
  # generating Betas
  foo <- inv(t(Z)%*%Z)
  Mu <- foo%*%t(Z)%*%y
  Sigma <- foo * sigma.2
  Betas <- mvrnorm(mu = Mu, Sigma = Sigma)
  return(list('sigma'= sigma.2 , 'Betas'=Betas))
}


gamma.distribution <- function(Z,y){
  n <- dim(Z)[1]
  k <- dim(Z)[2]
  inv <- det(t(Z)%*%Z)^(-1/2)
  foo <- (t(y)%*%(diag(n)- Z%*%inv(t(Z)%*%Z)%*%t(Z))%*%y)^(-(n-k-1)/2)
  print(c('inv_f'=inv,'foo_f'=foo))
  return(exp(inv) + exp(foo))
}

sampling_gammas <- function(gamma.matrix, Z, y, X){
  k <- dim(gamma.matrix)[1]
  p <- dim(gamma.matrix)[2]
  for(j in 1:k){
    new.G <- gamma.matrix
    new.gamma <- mvrnorm(mu = gamma.matrix[j,], Sigma = 0.05^2*diag(p))
    new.G[j,] <- new.gamma
    Z.new <- compute_Z_matrix(gamma.matrix = new.G, X.matrix = X)
    foo <- det(t(Z.new)%*%Z.new)
    #cat(c('foo'=foo, '\n'))
    if(foo>0.0001 & all(new.gamma< 1e5)){
      prob <- min(exp(1), gamma.distribution(Z.new,y) - gamma.distribution(Z,y))
      print(prob)
      if(prob <=0 ) break
      accept <- rbinom(1,1,prob = log(prob))
      if( accept==1){
        Z <- Z.new
        gamma.matrix <- new.G
      }
    }
  }
  return(list('gammas' = gamma.matrix, 'Z'= Z))
}




#  try --------------------------------------------------------------------

# #dat <- read.csv('dati_computer.txt', sep = ',', header = F)
# dat <- read.csv('MLB2008.csv', sep=',')
# dat <- dat[,-c(1,2,which(colnames(dat)=="PA_PR"))]
# normalizing <- function(x) (x- min(x))/(max(x)- min(x))
# dat[,2:131] <- apply(dat[,2:131],2,normalizing)
# X <- dat[,2:131]
# y <- dat[,1]
# 
# f <- as.formula(paste("SALARY~" , paste(names(dat)[!names(dat)=='SALARY'], collapse = '+')))
# 
# net.1 <- neuralnet(f, data = dat, hidden = 80, linear.output = T, lifesign = 'full', threshold = 0.1)
# net.dat <- neuralnet(f ,data = dat, hidden = 60, linear.output = T,
#                      lifesign = 'full', algorithm = 'backprop', learningrate = 0.1, threshold = 0.01, stepmax = 500, likelihood = T)
# plot(net.dat)
# 
# y_hat <- compute(net.dat, X)
# plot(y_hat$net.result , dat$Y)
# 
# 
# 
# 
# # MCMC
# K = 30
# p = dim(X)[2] + 1
# X.matrix<- as.matrix(X)
# y.vector <- as.array(y)
# 
# r = 1
# r
# gamma.matrix <- matrix(runif(K*p, -r ,r),c(K,p)) # la scelta difficile è quella di trovare un adeguato stato iniziale
# Z <- compute_Z_matrix(gamma.matrix, X.matrix)
# det(t(Z)%*%Z)
# 
# gamma.matrix <- MH_sampling$gammas
# Z <- MH_sampling$Z
# 
# conjugate_sampling <- sampling_sigma_and_betas(Z,y.vector) 
# conjugate_sampling$Betas
# MH_sampling <- sampling_gammas(gamma.matrix = gamma.matrix, Z = Z,y = y.vector,X = X.matrix)
# all(gamma.matrix==MH_sampling$gammas)
