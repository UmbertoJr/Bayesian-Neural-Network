# data <- read.csv('MLB2008.csv', sep=',')
# data <- data[, -c(1,2)]
# 
# observed.data <- as.matrix(data)
# X <- list(t(observed.data[,2:132]))
# dput(X, 'X_data.txt')
# 
# Y <- list(observed.data[,1])
# dput(Y, 'Y_data.txt')

dput(list(t(X)), 'X_data2.txt')
  
dput(list(Y), 'Y_data2.txt')

# initialization

Gamma = net.1$weights[[1]][[1]][-1,]
Gamma.0 = net.1$weights[[1]][[1]][1,]

Beta = net.1$weights[[1]][[2]][-1]
Beta.0 = net.1$weights[[1]][[2]][1]

dput(list(Gamma = t(Gamma), Gamma.0=Gamma.0, Beta=Beta, Beta.0=Beta.0), 'init_data2.txt')


  
# Building hyperparameters for first model --------------------------------

alpha.0 <- 1
alpha.1 <- 1
alpha.2 <- 1
alpha.a <- 1
alpha.b <- 1
omega.0 <- 10
omega.a <- 10
omega.b <- 10
omega.y <- 10

data <- list('alpha'= alpha.2/2, 'beta'= alpha.2/2/omega.y, 'alpha.o'= alpha.0/2, 'beta.o'= alpha.0/2/omega.0,
             'alpha.b'= alpha.b/2, 'beta.b'= alpha.b/2/omega.b, 'alpha.a'= alpha.a/2, 'beta.a'= alpha.a/2/omega.a,
             'alpha.out.1'= alpha.1/2, 'alpha.in.1'= alpha.1/2)

dput(data, 'hyper_parameters_first_model.txt')

# init priors

tau <- 0.15
tau.out <- 0.15
tau.in <- 0.15
tau.bo <- 0.15
tau.bin <- 0.15
tau.ao <- 0.15
tau.a <- rep(0.15, 20)
tau.out.vector <- rep(0.15, 20)
tau.in.vector <- rep(0.15, 5)

init_first_model <- list(tau = tau, tau.out = tau.out, tau.in=tau.in, tau.bo=tau.bo, tau.bin=tau.bin, tau.ao=tau.ao, tau.a=tau.a ,
                         tau.out.vector = tau.out.vector ,tau.in.vector= tau.in.vector)

dput(init_first_model , 'init_first_model.txt')





# Building hyperparameters for second model -------------------------------

R=dim(X)[2] + 1

a.beta = 0
A.beta = 1

a.gamma = rep(0, R)
A.gamma = diag(R)

alpha.beta = 0.5
beta.beta = 0.5

df = R + 1
scale = diag(R)

alpha.tau = 0.5
beta.tau = 0.05

dput(list(a.beta=a.beta, A.beta=A.beta, a.gamma=a.gamma, A.gamma=t(A.gamma), alpha.beta=alpha.beta, beta.beta=beta.beta,
          df=df, scale=t(scale), alpha.tau=alpha.tau, beta.tau=beta.tau, R = R),
          'hyper_parameters_second_model.txt')
#obs
X.hat <- cbind(rep(1,1000),X)
dput(list(t(X.hat)), 'X_data2_for_muller.txt')

# Init second model

list(Gamma=t(net.1$weights[[1]][[1]][,]),
     Beta.0 = net.1$weights[[1]][[2]][1],Beta = net.1$weights[[1]][[2]][-1])
dput(list(Gamma=t(net.1$weights[[1]][[1]][,]),
          Beta.0 = net.1$weights[[1]][[2]][1],Beta = net.1$weights[[1]][[2]][-1]),
     'init_data2_second_model.txt')

