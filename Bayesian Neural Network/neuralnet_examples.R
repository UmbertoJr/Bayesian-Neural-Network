library(neuralnet)


# examples  ---------------------------------------------------------------

# sum AND + OR
AND <- c(rep(0,7),1)
OR <- c(0, rep(1,7))
binary.data <- data.frame(expand.grid(c(0,1),c(0,1),c(0,1)), AND, OR)

net <- neuralnet(AND+OR~ Var1 + Var2 + Var3, binary.data, hidden = 0, rep = 10, err.fct = 'ce', linear.output = FALSE )
plot(net)
print(net)


# XOR
XOR <- c(0,1,1,0)
xor.data <- data.frame(expand.grid(c(0,1), c(0,1)),XOR)
net.xor <- neuralnet(XOR ~ Var1 + Var2, xor.data, hidden = 2, rep=5)
print(net.xor)
plot(net.xor, rep='best')


# infert data
data("infert", package = 'datasets')
net.infert <- neuralnet(case ~ parity + induced + spontaneous, infert, err.fct = 'ce', linear.output = FALSE, likelihood = T)
print(net.infert)
plot(net.infert)
gwplot(net.infert, selected.covariate = 'parity')
gwplot(net.infert, selected.covariate = 'induced')
gwplot(net.infert, selected.covariate = 'spontaneous')
confidence.interval(net.infert)


# Sqrt data
Var1 <- runif(50,0,100)
sqrt.data <- data.frame(Var1, Sqrt=sqrt(Var1))

net.sqrt <- neuralnet( Sqrt~ Var1,sqrt.data, hidden = 10, threshold = 0.01  )
plot(net.sqrt)
compute(net.sqrt, (1:10)^2)$net.result


# SUM data
Var1 <- rpois(100, 0.5)
Var2 <- rbinom(100, 2, 0.6)
Var3 <- rbinom(100,1, 0.5)
SUM <- as.integer(as.integer(abs(Var1+Var2+Var3+ rnorm(100))))
sum.data <- data.frame(Var1+Var2+Var3, SUM)
net.sum <- neuralnet(SUM ~Var1+Var2+Var3, sum.data, hidden = 1, act.fct = 'tanh' )
print(net.sum)
prediction(net.sum)
main <- glm(SUM ~Var1+Var2+Var3, sum.data, family = poisson())
full <- glm(SUM ~Var1*Var2*Var3, sum.data, family = poisson())
prediction(net.sum, list.glm = list(main=main, full=full))
