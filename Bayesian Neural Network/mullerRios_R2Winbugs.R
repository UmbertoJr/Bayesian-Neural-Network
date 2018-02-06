require(R2WinBUGS, quietly = T)
require(neuralnet, quietly = T)


model.file <- 'MullerRiosModel.txt'

file.show(model.file)

#load data
load('data_simmulated_for_neal.RData')

# DAta and hyper-paramaters

X=X
y = Y
N = dim(X)[1]
R=dim(X)[2] + 1

K=15

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

data <- list('X','y','N','R','K','a.beta','A.beta',
             'a.gamma','A.gamma','alpha.beta','beta.beta','df','scale',
             'alpha.tau','beta.tau')


# initialiazion

f <- as.formula(paste("Y~" , paste(names(dat)[!names(dat)=='Y'], collapse = '+')))
net.1 <- neuralnet(f, data = dat, hidden = K , linear.output = T,
                   lifesign = 'full', threshold = 0.1)


inits <- function(){
  list(Gamma=t(net.1$weights[[1]][[1]][-1,]),Gamma.0 = net.1$weights[[1]][[1]][1,],
       Beta.0 = net.1$weights[[1]][[2]][1],Beta = net.1$weights[[1]][[2]][-1])
}

parameters <- c('Gamma','Gamma.0','Beta.0','Beta','Y.rep','p.log','deviance')



# To run:

MullerRios_sim <-  bugs(data, inits, parameters, model.file,
                       n.chains = 3, n.iter = 5000,
                       bugs.directory = 'C:/Users/Umbertojunior/Desktop/applicazioni/applicazioni usate/statistica/winbugs14_unrestricted/WinBUGS14/',
                       debug = TRUE)
