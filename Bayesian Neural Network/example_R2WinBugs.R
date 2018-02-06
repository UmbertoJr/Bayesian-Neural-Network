require(R2WinBUGS, quietly = T)

## Same "schoolsmodel" that is used in the examples in ?bugs:
schoolsmodel <- function(){
  for (j in 1:J){
    y[j] ~ dnorm (theta[j], tau.y[j])
    theta[j] ~ dnorm (mu.theta, tau.theta)
    tau.y[j] <- pow(sigma.y[j], -2)
  }
  mu.theta ~ dnorm (0.0, 1.0E-6)
  tau.theta <- pow(sigma.theta, -2)
  sigma.theta ~ dunif (0, 1000)
}

filename <- file.path(tempdir(), "schoolsmodel.bug")

## write model file:
write.model(schoolsmodel, filename)
## and let's take a look:
file.show(filename)

# An example model file is given in:
model.file <- system.file(package = "R2WinBUGS", "model", "schools.txt")
# Let's take a look:
file.show(model.file)
# Some example data (see ?schools for details):
data(schools)
schools
J <- nrow(schools)

y <- schools$estimate
sigma.y <- schools$sd
data <- list ("J", "y", "sigma.y")
inits <- function(){
  list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100),
       sigma.theta = runif(1, 0, 100))
}
## or alternatively something like:
# inits <- list(
#   list(theta = rnorm(J, 0, 90), mu.theta = rnorm(1, 0, 90),
#        sigma.theta = runif(1, 0, 90)),
#   list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100),
#        sigma.theta = runif(1, 0, 100))
#   list(theta = rnorm(J, 0, 110), mu.theta = rnorm(1, 0, 110),
#        sigma.theta = runif(1, 0, 110)))

parameters <- c("theta", "mu.theta", "sigma.theta")
## Not run: 
## both write access in the working directory and package BRugs required:
schools.sim <- bugs(data, inits, parameters, model.file,
                    n.chains = 3, n.iter = 5000,
                    bugs.directory = 'C:/Users/Umbertojunior/Desktop/applicazioni/applicazioni usate/statistica/winbugs14_unrestricted/WinBUGS14/')
print(schools.sim)
plot(schools.sim)
