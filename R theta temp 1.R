rm(list=ls())
library(boot)

#===============
# simulate data
#===============
nl <- 100
nc <- 2
ns <- 5
nmat <- 10

niter <- 10000
nburn <- round(niter / 2)

gamma <- .5

vmat <- matrix(rbeta(nl*nc, 1, gamma), nl, nc)
vmat[,nc] <- 1

vmat2theta <- function(vmat) {
  nl <- dim(vmat)[1]
  nc <- dim(vmat)[2]
  theta <- vmat
  for (i in 2:nc) {
    for (j in 1:nl) {
      theta[j,i] <- vmat[j,i] * cumprod(1-vmat[j,1:(i-1)])[i-1]
    } # j
  } # i
  return(theta)
} # function

theta <- vmat2theta(vmat)
colMeans(vmat)
colMeans(theta)

phi <- matrix(rgamma(nc*ns, 5, 1), nc, ns)
for (i in 1:nc) {
  phi[i,] <- phi[i,] / sum(phi[i,])
}

pi <- theta %*% phi

y <- matrix(rbinom(nl*ns, nmat, pi), nl, ns)
head(dbinom(y, nmat, pi), 5)

param <- list()
param$phi <- phi
param$theta <- theta
param$vmat <- vmat

#==================
# define functions
#==================
tnorm <- function(n, lo, hi, mu, sig) {
  # generates truncated normal variates based on cumulative normal distribution normal truncated lo and hi

  if (length(lo) == 1 & length(mu) > 1) 
    lo <- rep(lo, length(mu))
  if (length(hi) == 1 & length(mu) > 1) 
    hi <- rep(hi, length(mu))
    
  q1 <- pnorm(lo, mu, sig)  #cumulative distribution
  q2 <- pnorm(hi, mu, sig)  #cumulative distribution
    
  z <- runif(n, q1, q2)
  z <- qnorm(z, mu, sig)
  z[z == -Inf] <- lo[z == -Inf]
  z[z == Inf] <- hi[z == Inf]
  z
}

fix.MH <- function(lo, hi, old1, new1, jump) {
  jold <- pnorm(hi, mean = old1, sd = jump) - pnorm(lo, mean = old1, sd = jump)
  jnew <- pnorm(hi, mean = new1, sd = jump) - pnorm(lo, mean = new1, sd = jump)
  log(jold) - log(jnew)  #add this to pnew
}

fix.probs = function(probs) {
    cond = probs < 1e-05
    probs[cond] = 1e-05
    cond = probs > 0.99999
    probs[cond] = 0.99999
    probs
}

get.logl = function(theta, phi, y, nmat) {
    prob = fix.probs(theta %*% phi)
    dbinom(y, size = nmat, prob = prob, log = T)
}

acceptMH <- function(p0, p1, x0, x1, BLOCK) {
  # accept for M, M-H if BLOCK, then accept as a block, otherwise, accept individually
    
  nz <- length(x0)  #no. to accept
  if (BLOCK) 
    nz <- 1
    
  a <- exp(p1 - p0)  #acceptance PR
  z <- runif(nz, 0, 1)
  keep <- which(z < a)
    
  if (BLOCK & length(keep) > 0) 
    x0 <- x1
  if (!BLOCK) 
    x0[keep] <- x1[keep]
  accept <- length(keep)
    
  list(x = x0, accept = accept)
}

update.theta <- function(param, jump, nl, nc, y, nmat, gamma) {
    v.old <- param$vmat
    v.tmp <- tnorm(nl * (nc - 1), lo = 0, hi = 1, mu = v.old[, -nc], sig = jump)
    novos <- cbind(matrix(v.tmp, nl, nc - 1), 1)
    adj <- matrix(fix.MH(lo = 0, hi = 1, old1 = v.old, new1 = novos, jump = jump), nl, nc)
    
    prior.old <- matrix(dbeta(v.old, 1, gamma, log = T), nl, nc)
    prior.new <- matrix(dbeta(novos, 1, gamma, log = T), nl, nc)
    
    for (j in 1:(nc-1)) {
        # last column has to be 1
        v.new <- v.old
        v.new[,j] <- novos[,j]
        
        theta.old <- vmat2theta(vmat = v.old)
        theta.new <- vmat2theta(vmat = v.new)
        
        pold <- get.logl(theta = theta.old, phi = param$phi, y = y, nmat = nmat)
        pnew <- get.logl(theta = theta.new, phi = param$phi, y = y, nmat = nmat)
        p1.old <- rowSums(pold)
        p1.new <- rowSums(pnew)
        
        k <- acceptMH(p0=p1.old+prior.old[,j], p1=p1.new+prior.new[,j]+adj[,j], x0=v.old[,j], x1=v.new[,j], F)
        v.old[,j] <- k$x
    }
    theta <- vmat2theta(vmat=v.old)
    list(theta=theta, vmat=v.old)
}

#==================
# test
#==================
theta.post <- array(, dim=c(nl, nc, niter))
theta.post[,,1] <- param$theta
for (i in 2:niter) {
  theta.up <- update.theta(param, jump=.1, nl, nc, y, nmat, gamma=gamma)
  param$theta <- theta.up$theta
  param$vmat <- theta.up$vmat
  theta.post[,,i] <- param$theta
}

theta.est <- apply(theta.post[,,nburn:niter], 1:2, median)

par(mfrow=c(3,2))
plot(theta.post[1,1,], type='l'); abline(h=theta[1,1], col=2)
plot(theta.post[1,2,], type='l'); abline(h=theta[1,2], col=2)
plot(theta.post[2,1,], type='l'); abline(h=theta[2,1], col=2)
plot(theta.post[2,2,], type='l'); abline(h=theta[2,2], col=2)

plot(theta.est ~ theta, xlim=range(c(theta,theta.est)), ylim=range(c(theta,theta.est))); abline(0, 1, col=2)

boxplot(theta.est)


