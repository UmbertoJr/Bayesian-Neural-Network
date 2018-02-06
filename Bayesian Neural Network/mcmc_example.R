library(mcmc)
data("logit")

out <- glm(y~ x1 + x2 + x3+ x4, data = logit, family = binomial(), x=T)
summary(out)



# Bayesian setup ----------------------------------------------------------

x <- out$x
y <- out$y

lupost <- function(beta, x, y) {
  eta <- as.numeric(x %*% beta)
  logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
  logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    logl <- sum(logp[y == 1]) + sum(logq[y == 0])
  return(logl - sum(beta^2) / 8)
  }


# beggining MCMC ----------------------------------------------------------

set.seed(42)
beta.int <- as.numeric(coefficients(out))
out <- metrop(lupost, beta.int, 1e3, x=x, y=y)
names(out)

out$accept

# find the right scale
out <- metrop(out, scale = 0.1, x=x, y=y)
out$accept

out <- metrop(out, scale = 0.3, x=x, y=y)
out$accept

out <- metrop(out, scale = 0.5, x=x, y=y)
out$accept

out <- metrop(out, scale = 0.4, x=x, y=y)
out$accept


# Diagnostic --------------------------------------------------------------


out <- metrop(out, nbatch = 1e4, x=x, y=y)
out$accept
out$time

plot(ts(out$batch))
acf(out$batch)



# Monte Carlo Estimates and Standard Errors -------------------------------

out <- metrop(out, nbatch = 1e2, blen = 100, outfun = function(z, ...)c(z,z^2), x=x, y=y) # bacthes of 100
out$accept
out$time

# Simple Means

foo <- apply(out$batch, 2, mean)
mu <- foo[1:5]
sigmasq <- foo[6:10] - mu^2
mu
sigmasq

mu.mcse <- apply(out$batch[, 1:5], 2 , sd) / sqrt(out$nbatch)
mu.mcse

# Variance MCSE

u <- out$batch[, 1:5]
v <- out$batch[,6:10]
ubar <- apply(u, 2, mean)
vbar <- apply(v, 2, mean)
deltau <- sweep(u, 2, ubar)
deltav <- sweep(v, 2, vbar)
foo <- sweep(deltau, 2, ubar,'*')
sigmasq.mcse <- sqrt(apply((deltau - 2*foo)^2, 2,mean)/ out$nbatch)
sigmasq.mcse

sigma <- sqrt(sigmasq)
sigma.mcse <- sigmasq.mcse / (2*sigma)
sigma
sigma.mcse




# New Variance Estimator --------------------------------------------------

n <- 2e4
rho <- 0.99
x <- arima.sim(model = list(ar=rho), n=n)

out <- initseq(x)
names(out)

plot(seq(along = out$Gamma.pos)-1 , out$Gamma.pos, xlab = 'k', ylab = expression(Gamma[k]),type = 'l')
lines(seq(along= out$Gamma.dec)- 1, out$Gamma.dec, lty= 'dotted')
lines(seq(along= out$Gamma.con)- 1, out$Gamma.con, lty= 'dashed')

out$var.con

(1 - rho)/ (1 - rho) * 1 /(1 - rho)^2

blen <- 5
x.batch <- apply( matrix(x , nrow = blen), 2 , mean)
bout <- initseq(x.batch)

plot(seq(along= bout$Gamma.con) - 1, bout$Gamma.con , xlab = 'k', ylab = expression(Gamma[k]),type = 'l' )








