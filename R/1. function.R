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
} # tnorm

fix.MH <- function(lo, hi, old1, new1, jump) {
  jold <- pnorm(hi, mean = old1, sd = jump) - pnorm(lo, mean = old1, sd = jump)
  jnew <- pnorm(hi, mean = new1, sd = jump) - pnorm(lo, mean = new1, sd = jump)
  log(jold) - log(jnew)  #add this to pnew
} # fix.MH

get.logl <- function(theta, phi, pobs.mat, y, nmat) {
  prob <- theta %*% phi * pobs.mat
  cond <- prob < 1e-05
  prob[cond] <- 1e-05
  cond <- prob > 0.99999
  prob[cond] <- 0.99999
  dbinom(y, size = nmat, prob = prob, log = T)
} # get.logl

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
} # acceptMH

update.theta <- function(param, jump, nl, nc, ns, y, x, nmat, id.perf) {
  phi <- param$phi
  beta <- param$beta
  mu <- x %*% beta
  pobs <- param$pobs
  pobs.mat <- matrix(pobs, nl, ns, byrow=T)
  pobs.mat[id.perf,] <- 1

  vmat.ori <- vmat.old <- param$vmat
  vmat.tmp <- exp(rnorm(nl*(nc-1), mean=log(vmat.old[,-nc]), sd=jump[,-nc]))
  proposed <- cbind(matrix(vmat.tmp, nl, nc-1), 1)

  for (j in 1:(nc-1)) {
    # last column has to be 1
    vmat.new <- vmat.old
    vmat.new[,j] <- proposed[,j]

    theta.old <- vmat.old / rowSums(vmat.old)
    theta.new <- vmat.new / rowSums(vmat.new)

    prob.old <- get.logl(theta=theta.old, phi=phi, pobs.mat=pobs.mat, y=y, nmat=nmat)
    prob.new <- get.logl(theta=theta.new, phi=phi, pobs.mat=pobs.mat, y=y, nmat=nmat)

    pold <- rowSums(prob.old) + dnorm(log(vmat.old[,j]), mu[,j], sigma, log = T)
    pnew <- rowSums(prob.new) + dnorm(log(vmat.new[,j]), mu[,j], sigma, log = T)

    k <- acceptMH(p0=pold, p1=pnew, x0=vmat.old[,j], x1=vmat.new[,j], BLOCK=F)
    vmat.old[,j] <- k$x
  }
  vmat <- vmat.old
  theta <- vmat / rowSums(vmat)
  list(theta=theta, vmat=vmat, accept=vmat.ori!=vmat.old)
} # update.theta

update.phi <- function(param, jump, nl, nc, ns, y, nmat, id.perf, a.phi, b.phi) {
  theta <- param$theta
  pobs <- param$pobs
  pobs.mat <- matrix(pobs, nl, ns, byrow=T)
  pobs.mat[id.perf,] <- 1

  phi.ori <- phi.old <- param$phi
  proposed <- matrix(tnorm(nc*ns, lo=0, hi=1, mu=phi.old, sig=jump), nc, ns)
  adj <- fix.MH(lo=0, hi=1, old1=phi.old, new1=proposed, jump=jump)

  for (j in 1:nc) {
    phi.new <- phi.old
    phi.new[j,] <- proposed[j,]

    prob.old <- get.logl(theta=theta, phi=phi.old, pobs.mat=pobs.mat, y=y, nmat=nmat)
    prob.new <- get.logl(theta=theta, phi=phi.new, pobs.mat=pobs.mat, y=y, nmat=nmat)

    pold <- colSums(prob.old) + dbeta(phi.old[j,], a.phi, b.phi, log=T)
    pnew <- colSums(prob.new) + dbeta(phi.new[j,], a.phi, b.phi, log=T)

    k <- acceptMH(p0=pold, p1=pnew + adj[j,], x0=phi.old[j,], x1=phi.new[j,], BLOCK=F)
    phi.old[j,] <- k$x
  }
  phi <- phi.old
  list(phi=phi, accept=phi.ori!=phi.old)
} # update.phi

update.lpmu <- function(param, jump, ns) {
  pobs <- param$pobs

  lpmu.old <- param$lpmu
  lpsigma.old <- param$lpsigma

  lpmu.new <- rnorm(1, lpmu.old, jump)
  lpsigma.new <- exp(rnorm(1, log(lpsigma.old), jump))

  pold <- sum(dnorm(logit(pobs), lpmu.old, lpsigma.old, log=T)) + 
          dnorm(lpmu.old, 0, 10, log=T) + dnorm(log(lpsigma.old), 0, 10, log=T)
  pnew <- sum(dnorm(logit(pobs), lpmu.new, lpsigma.new, log=T)) + 
          dnorm(lpmu.new, 0, 10, log=T) + dnorm(log(lpsigma.new), 0, 10, log=T)

  k <- acceptMH(p0=pold, p1=pnew, x0=c(lpmu.old,lpsigma.old), x1=c(lpmu.new,lpsigma.new), BLOCK=F)
  lpmu <- k$x[1]
  lpsigma <- k$x[2]
  list(lpmu=lpmu, lpsigma=lpsigma)
} # update.lpmu

update.pobs <- function(param, jump, nl, ns, id.perf) {
  theta <- param$theta
  phi <- param$phi
  lpmu <- param$lpmu
  lpsigma <- param$lpsigma

  pobs.old <- param$pobs
  pobs.new <- inv.logit(rnorm(ns, logit(pobs.old), jump))
  pobs.new <- ifelse(pobs.new > .99, .99, pobs.new)
  pobs.new <- ifelse(pobs.new < .01, .01, pobs.new)
  pobs.old.mat <- matrix(pobs.old, nl, ns, byrow=T)
  pobs.old.mat[id.perf,] <- 1
  pobs.new.mat <- matrix(pobs.new, nl, ns, byrow=T)
  pobs.new.mat[id.perf,] <- 1

  prob.old <- get.logl(theta=theta, phi=phi, pobs.mat=pobs.old.mat, y=y, nmat=nmat)
  prob.new <- get.logl(theta=theta, phi=phi, pobs.mat=pobs.new.mat, y=y, nmat=nmat)

  pold <- sum(prob.old) + sum(dnorm(logit(pobs.old), lpmu, lpsigma, log=T))
  pnew <- sum(prob.new) + sum(dnorm(logit(pobs.new), lpmu, lpsigma, log=T))

  k <- acceptMH(p0=pold, p1=pnew, x0=pobs.old, x1=pobs.new, BLOCK=F)
  pobs <- k$x
  list(pobs=pobs, accept=pobs!=pobs.old)
} # update.pobs

update.beta <- function(param, jump, nl, nc, x) {
  vmat <- param$vmat

  beta.old <- param$beta
  beta.new <- matrix(rnorm(length(beta.old), beta.old, jump), nx, nc)
  mu.old <- x %*% beta.old
  mu.new <- x %*% beta.new

  prob.old <- dnorm(log(vmat[,-nc]), mu.old[,-nc], sigma, log=T)
  prob.new <- dnorm(log(vmat[,-nc]), mu.new[,-nc], sigma, log=T)

  pold <- sum(prob.old) + sum(dnorm(beta.old[-nc], 0, 5, log=T))
  pnew <- sum(prob.new) + sum(dnorm(beta.new[-nc], 0, 5, log=T))

  k <- acceptMH(p0=pold, p1=pnew, x0=beta.old, x1=beta.new, BLOCK=F)
  beta <- k$x
  list(beta=beta, accept=beta!=beta.old)
} # update.beta

thetaphi.jumpTune <- function(accept, jump, ni, adapt=2000, low=.3, high=.8) {
  nstart <- ifelse(ni>=100, ni-99, 1)
  if (ni > adapt) {
    jump.new <- jump
  } else if (ni %% 10 > 0) {
    jump.new <- jump
  } else {
    accept.rate <- apply(accept[,,nstart:ni], 1:2, mean)
    jump.new <- jump
    jump.new[which(accept.rate < low)] <- 
      jump[which(accept.rate < low)] * rnorm(length(which(accept.rate < low)),.5,.01)
    jump.new[which(accept.rate > high)] <- 
      jump[which(accept.rate > high)] * rnorm(length(which(accept.rate > high)),2,.01)
  }
  jump <- jump.new
  list(jump=jump)
} # jumpTune

pobsbeta.jumpTune <- function(accept, jump, ni, adapt=2000, low=.3, high=.8) {
  nstart <- ifelse(ni>=100, ni-99, 1)
  if (ni > adapt) {
    jump.new <- jump
  } else if (ni %% 10 > 0) {
    jump.new <- jump
  } else {
    accept.rate <- mean(accept[nstart:ni], na.rm=T)
    if (accept.rate < low) {
      jump.new <- jump * rnorm(1,.5,.05)
    } else if (accept.rate > high) {
      jump.new <- jump * rnorm(1,2,.05)
    } else {
      jump.new <- jump
    }
  }
  jump <- jump.new
  list(jump=jump)
} # pobsbeta.jumpTune


