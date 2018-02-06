require(rjags)
require(coda)


# Normally distributed data with unknow mean and variance -----------------

set.seed(432104)
n <- 1000
x <- rnorm(n, 0,5)


model1.string <- "
  model {
    for(i in 1:N){
    x[i] ~ dnorm(mu, tau)
    }
  mu ~ dnorm(0, 0.0001)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
  }
"

model1.spec <- textConnection(model1.string)

jags <- jags.model(model1.spec, 
                   data = list('x'=x, 
                               'N'=n),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 1e3)

jags.samples(jags, 
             c('mu', 'tau'),
             1000)



# SAT test example --------------------------------------------------------

sigma     <- c(15,10,16,11, 9,11,10,18)
schoolobs <- c(28,8, -3, 7,-1, 1,18,12)


model.sat.text<-"
  model {
for(i in 1:N) {
schoolmean[i] ~ dnorm(mu,itau)
thes[i] <- 1/pow(sigma[i],2)
schoolobs[i] ~ dnorm(schoolmean[i],thes[i])
}

mu ~ dnorm(0,alpha)
alpha <- .01
itau   ~ dgamma(1e-3,pow(15,2)*1e-3)
tau <- pow(1/itau,1/2)
}
"

model.sat.spec<-textConnection(model.sat.text)

require(igraph)

gr<-graph.formula("N(0,0.01)"-+"mu",
                  "mu"-+"N(0,1/tau)", 
                  "N(0,1/tau)"-+"m1", 
                  "N(0,1/tau)"-+"m2", 
                  "N(0,1/tau)"-+"m8",
                  "m1"-+"N(0,1/simga21)", 
                  "m2"-+"N(0,1/simga22)",
                  "m8"-+"N(0,1/simga28)", 
                  "N(0,1/simga21)"-+"y1",
                  "N(0,1/simga22)"-+"y2",  
                  "N(0,1/simga28)"-+"y8")


lo<-data.frame(x=c(2,2,2,1,2,3,1,2,3,1,2,3),y=c(6,5,4,3,3,3,2,2,2,1,1,1))
plot(gr, 
     layout=layout.reingold.tilford(gr), 
     edge.arrow.size=.25
)


sat.jags <- jags.model(model.sat.spec,
                       data=list('sigma'=sigma,
                                 'schoolobs'=schoolobs,
                                 'N'=length(schoolobs)
                       ),
                       n.adapt = 1000)

samps.jags <- jags.samples(sat.jags,
                           c('mu','tau'),
                           n.iter=10000,
                           thin=10
)

samps.jags


samps.coda <- coda.samples(sat.jags,
                           c('mu','tau', 'schoolmean'),
                           n.iter=10000,
                           thin=10
)

head(samps.coda)

summary(samps.coda)


plot(samps.coda[[1]][,c("mu")])

class(samps.coda[[1]])



# Linear regression example -----------------------------------------------

N <- 1000
x <- 1:N
epsilon <- rnorm(N,0,1)
y <- x + epsilon

write.table(data.frame(X = x, Y = y, Epsilon = epsilon),
            file = 'example2.data',
            row.names = FALSE,
            col.names = TRUE)


jags <- jags.model('example2.bug',
                   data = list('x' = x,
                               'y' = y,
                               'N' = N),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 1000)

jags.samples(jags,
             c('a', 'b'),
             1000)



# Logistical regression example -------------------------------------------

N <- 1000
x <- 1:N
z <- 0.01 * x - 5
y <- sapply(1 / (1 + exp(-z)), function(p) {rbinom(1, 1, p)})

write.table(data.frame(X = x, Z = z, Y = y),
            file = 'example3.data',
            row.names = FALSE,
            col.names = TRUE)

jags <- jags.model('example3.bug',
                   data = list('x' = x,
                               'y' = y,
                               'N' = N),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 1000)

jags.samples(jags,
             c('a', 'b'),
             1000)
