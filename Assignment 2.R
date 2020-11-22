## Question 2

load("dataex2.Rdata")
require(maxLik)

log_lik <- function(mu, data) {
  x <- data[,1]; r <- data[,2]
  sigma <- 1.5
  sum(r*dnorm(x, mu, sigma, log=TRUE) + (1-r)*pnorm(x, mu, sigma, log.p=TRUE))
}

mle <- maxLik(logLik = log_lik, data = dataex2, start = c(mu = 1))
summary(mle)


## Question 4
load("dataex4.Rdata")
ind <- which(is.na(dataex4$Y) == FALSE)
y.obs <- dataex4$Y[ind]
x.m <- dataex4$X[ind]
x.n_m <- dataex4$X[-ind]
params.old <- c(0,0)

p <- function(params, x){
  beta0 <- params[1]; beta1 <- params[2]
  exp(beta0 + x*beta1)/(1 + exp(beta0 + x*beta1))
}

Q <- function(params){
  beta0 <- params[1]; beta1 <- params[2]
  -(sum(y.obs*(beta0 + x.m*beta1) - log(1+exp(beta0 +x.m*beta1))) + sum(p(params.old, x.n_m)*(beta0 + x.n_m*beta1)))
}

em.optim <- function(params0, eps) {
  diff <- 1
  params <- params0
  while(diff > eps) {
    params.old <- params
    params <- optim(params.old, Q)$par
    diff <- sum(abs(params- params.old))
  }
  return(params)
}

res <- em.optim(c(0,0), 1e-8)
beta0 <- res[1]; beta1 <- res[2]
beta0; beta1

## Question 5
em.mixture <- function(y, theta0, eps){
  n <- length(y)
  theta <- theta0
  p <- theta[1]
  mu <- theta[2]; sigma <- theta[3]
  lambda <- theta[4]
  diff <- 1
  while(diff > eps){
    theta.old <- theta
    #E-step
    ptilde1 <- p*dlnorm(y, mean = mu, sd = sigma)
    ptilde2 <- (1 - p)*dexp(y, rate = lambda)
    ptilde <- ptilde1/(ptilde1 + ptilde2)
    #M-step
    p <- mean(ptilde)
    mu <- sum(log(y)*ptilde)/sum(ptilde)
    sigma <- sqrt(sum(((log(y) - mu)^2)*ptilde)/sum(ptilde))
    lambda <- sum(1- ptilde)/sum(y*(1-ptilde))
    theta <- c(p, mu, sigma, lambda)
    diff <- sum(abs(theta - theta.old))
  }
  return(theta)
}

load("dataex5.Rdata")
hist(dataex5, main = "Histogram of dataex5",
     xlab = "dataex5",
     ylab = "", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.4)

theta0 <- c(0.1, 1, 0.5, 2)
res <- em.mixture(y = dataex5, theta0, 0.00001)
p <- res[1]
mu <- res[2]
sigma <- res[3]
lambda <- res[4]

p; mu; sigma^2; lambda

hist(dataex5, main = "Histogram of dataex5",
     xlab = "dataex5",
     ylab = "Density",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.4,
     freq = F, ylim=c(0,0.18))
curve(p*dlnorm(x, mean = mu, sd = sigma) + (1 - p)*dexp(x, rate = lambda),
      add = TRUE, lwd = 3, col = "blue")
