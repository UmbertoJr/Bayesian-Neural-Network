require(rjags)
require(coda)
require(neuralnet)

# model neal for BNN ------------------------------------------------------


#Dataset creation
N=1000     # number of observations
P=5       # number of covariates
K = 20    # number of neurons in the hidden layer

# X <- matrix(runif(N*P), c(N,P))
# Y <- (X[,1]^2 + log(X[,2]) + 3^X[,3]+ X[,4]*X[,5])# + c(rnorm(100))
# dat <- data.frame(cbind(X,Y))
#save(file = 'data_simmulated_for_neal.RData', dat, X,Y)

# load data
#load(file = 'data_simmulate.RData')
load('data_simmulated_for_neal.RData')
# fitting normal NN
#K=18
K=15
f <- as.formula(paste("Y~" , paste(names(dat)[!names(dat)=='Y'], collapse = '+')))
net.1 <- neuralnet(f, data = dat, hidden = K , linear.output = T,
                   lifesign = 'full', threshold = 0.1)

save(net.1, file = 'network_for_K_15_1000.RData')
Z <- compute_Z_matrix(net.1$weights[[1]][[1]],X)

Y.hat <- Z%*%net.1$weights[[1]][[2]]
plot(net.1$net.result[[1]],-32 -Y.hat)

plot(net.1$net.result[[1]], Y)
sum((net.1$net.result[[1]]- Y)^2)/42


# Jags model compile
jags3 <- jags.model('neal2.bug', data = list('X'=X, 'y'=Y, 'K'=K, 'R'=5, 'N'=N,
                                           'alpha'= alpha.2/2, 'beta'= alpha.2/2/omega.y,
                                           'alpha.o'= alpha.0/2, 'beta.o'= alpha.0/2/omega.0,
                                           'alpha.b'= alpha.b/2, 'beta.b'= alpha.b/2/omega.b,
                                             'alpha.a'= alpha.a/2, 'beta.a'= alpha.a/2/omega.a,
                                           'alpha.out.1'= alpha.1/2, 'alpha.in.1'= alpha.1/2),
                    inits = list('Gamma.0' = net.1$weights[[1]][[1]][1,], 'Gamma'=t(net.1$weights[[1]][[1]][-1,]),
                                 'Beta'= net.1$weights[[1]][[2]][-1], 'Beta.0'= net.1$weights[[1]][[2]][1]),
                   n.adapt = 1000)


update(jags3, 100)

sim_neal_huge <- jags.samples(jags3,variable.names = c('Gamma', 'Beta',
                                                 'Gamma.0','Beta.0',
                                                 'tau', 'Y.rep','p.log'),100000)

save(sim_neal_huge, file = 'neal_chain_huge_chain.RData')

load('neal_chain2.RData')
load('neal_chain3_for_15nuerons.RData')
load('neal_chain4_for_15nuerons.RData')
load('neal_chain_huge_chain.RData')

# visualizing distribution ------------------------------------------------

sim <- sim_neal_huge

for(k in 1:K){
  plot_and_find_the_mode(sim$Gamma.0,k, obj='Gamma.0')
  Sys.sleep(1.5)
  for(p in 1:P){
    plot_and_find_the_mode(sim$Gamma,k,p, obj='Gamma')
    Sys.sleep(1.5)
  }
}

plot_and_find_the_mode(sim_neal2$Beta.0, obj = 'Beta.0')
for(k in 1:K){
  plot_and_find_the_mode(sim_neal2$Beta,k, obj = 'Beta')
  Sys.sleep(1.5)
}

for(i in 1:10){
  plot(density(sim_neal2$Y.rep[i,,1]))
  points(Y[i],0,col='blue')
  p <- mean(sim_neal2$Y.rep[i,,1]< Y[i])
  p_val <- min(p, 1-p)
  print(p_val)
  Sys.sleep(0.7)
}

# Traceplots

ts.plot(sim$Beta[1,,])
ts.plot(sim$Beta.0[1,,])
ts.plot(sim$Beta[2,,])
ts.plot(sim$Gamma[1,1,,])

# running mean

sim1 <- sim$Beta[1,,]
sim1 <- sim$Beta.0[1,,]
sim1 <- sim$Beta[2,,]
sim1 <- sim$Gamma[1,1,,]

ts.plot(cumsum(sim1)/seq_along(sim1))

# Make predictions

y.hat.mode.2 <- prediction_for_neal_BNN(X,sim = sim)
y.hat.mean.2 <- prediction_for_neal_BNN_with_mean(X,sim)

plot(Y,y.hat,xlim = c(-4,5),ylim = c(0,4))
plot(Y,y.hat.mean.2)
plot(Y,y.hat.mode.2)

# # analisys of goodness-of-fit
# Y.rep <- replicated_from_posterior_predictive(X, sim_neal2)
# save(Y.rep, file = 'Y_replicated_from_neal_priors.RData')

# function used -----------------------------------------------------------

plot_and_find_the_mode <- function(sim, i =0, j=0, obj=''){
  if(i == 0){
    den <- density(sim[,,1])
  }else if(j==0){
    den <- density(sim[i,,1])
  }else{
    den <- density(sim[i,j,,1])
  }
  massimo <- den$x[which.max(den$y)]
  if(obj=='Beta'){
    lim <- net.1$weights[[1]][[2]][i,1]
    plot(den, xlim=c(-abs(lim)-2,+ abs(lim)+2), main = paste('Beta',i) )
    segments(massimo,0,y1=max(den$y), col='grey')
    points(lim ,0, col='blue')
  }
  if(obj=='Beta.0'){
    lim <- net.1$weights[[1]][[2]][(K+1),1]
    plot(den, xlim=c(-abs(lim)-2,+ abs(lim)+2), main ='Beta0')
    segments(massimo,0,y1=max(den$y), col='grey')
    points(lim,0, col='blue')
  }
  if(obj=='Gamma'){
    lim <- net.1$weights[[1]][[1]][j,i]
    plot(den, xlim=c(-abs(lim)-2,+ abs(lim)+2),main = paste(paste('Gamma- neuron: ',i),
                                                           paste(' -covariate: ', j)) )
    segments(massimo,0,y1=max(den$y), col='grey')
    points(lim,0, col='blue')
  }
  if(obj=='Gamma.0'){
    lim <- net.1$weights[[1]][[1]][6,i]
    plot(den, xlim=c(-abs(lim)-2,+ abs(lim)+2), main = paste('Gamma.0 Bias - neuron: ',i))
    segments(massimo,0,y1=max(den$y), col='grey')
    points(lim,0, col='blue')
  }
  return(massimo)
}

find_maximum <- function(sim, i = 0, j=0){
  if(i == 0){
    den <- density(sim[,,1])
  }else if(j==0){
    den <- density(sim[i,,1])
  }else{
    den <- density(sim[i,j,,1])
  }
  massimo <- den$x[which.max(den$y)]
  return(massimo)
}

prediction_for_neal_BNN <- function(X,sim){
  bias.beta <- find_maximum(sim_neal$Beta.0)
  betas <- rep(NA,20)
  for(n in 1:20){
    betas[n]<- find_maximum(sim_neal$Beta,i=n)
  }
  bias.gammas <- rep(NA,20)
  gammas <- matrix(rep(NA,20*5),c(20,5))
  for(n in 1:20){
    bias.gammas[n] <- find_maximum(sim_neal$Gamma.0,i=n)
    for(p in 1:5){
      gammas[n,p]<- find_maximum(sim_neal$Gamma,i = n,j = p)
    }
  }
  N <- dim(X)[1]
  y.hat <- rep(NA,N)
  Z <- matrix(rep(NA,N*20),c(N,20))
  Z.np <- matrix(rep(NA,N*20),c(N,20))
  for(i in 1:N){
    for(j in 1:20){
      Z[ i, j] <- bias.gammas[ j ] +  X[i,]%*%gammas[j,]
      #print(Z[ i, j])
      Z.np[i,j] <-  (1) / (1+ exp(- Z[i,j]))
    }
    y.hat[i] <- bias.beta +Z.np[i,]%*%betas
  }
  return(y.hat)
}

prediction_for_neal_BNN_with_mean <- function(X,sim){
  bias.beta <- mean(sim_neal$Beta.0)
  betas <- rep(NA,20)
  for(n in 1:20){
    betas[n]<- mean(sim_neal$Beta[n,,1])
  }
  bias.gammas <- rep(NA,20)
  gammas <- matrix(rep(NA,20*5),c(20,5))
  for(n in 1:20){
    bias.gammas[n] <- mean(sim_neal$Gamma.0[n,,1])
    for(p in 1:5){
      gammas[n,p]<- mean(sim_neal$Gamma[n,p,,1])
    }
  }
  N <- dim(X)[1]
  y.hat <- rep(NA,N)
  Z <- matrix(rep(NA,N*20),c(N,20))
  Z.np <- matrix(rep(NA,N*20),c(N,20))
  for(i in 1:N){
    for(j in 1:20){
      Z[ i, j] <- bias.gammas[ j ] +  X[i,]%*%gammas[j,]
      Z.np[i,j] <-  (1) / (1+ exp(- Z[i,j]))
    }
    y.hat[i] <- bias.beta +Z.np[i,]%*%betas
  }
  return(y.hat)
}

replicated_from_posterior_predictive <- function(X, chain){
 L <- 10000
 d <- dim(X)[1]
 y.rep <- matrix(rep(NA,L*d),nrow = L, ncol = d)
 for(i in 1:L){
   betas <- c(chain$Beta.0[1,i,1], chain$Beta[,i,1])
   gamma <- rbind((chain$Gamma.0[,i,1]),t(chain$Gamma[,,i,1]))
   Z <- compute_Z_matrix(gamma,X)
   mu_y <- Z%*%betas
   sigma <- 1/sqrt(chain$tau[1,i,1])* diag(d)
   y.rep[i] <- mvrnorm(mu = mu_y, Sigma = sigma)
 }
 return(y.rep)
}


