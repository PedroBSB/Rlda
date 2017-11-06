rm(list=ls())
library(boot)

path <- 'c:/A. UFL/1. regression/1. simulation/'

#============================
# define regression function
#============================
rlda.regression <- function(y=y, x=x, nmat=nmat, id.binomial=id.binomial, nc=nc, 
                            niter=niter, nburn=nburn, nadapt=nadapt, 
                            vmat.init=vmat.init, phi.init=phi.init, beta.init=beta.init, 
                            logit_pobs_mu.init=logit_pobs_mu.init, 
                            logit_pobs_sigma.init=logit_pobs_sigma.init, 
                            pobs.init=pobs.init) {
#======================
# call other functions
#======================
source(paste(c(path, 'code/1. function.R'), collapse=''))

nl <- dim(y)[1] # number of locations
ns <- dim(y)[2] # number of species
nx <- dim(x)[2] # number of covariates

#===========
# run model
#===========
theta.jump <- matrix(.05, nrow=nl, ncol=nc)
phi.jump <- matrix(.05, nrow=nc, ncol=ns)
pobs.jump <- .05
beta.jump <- .05

param <- list()

param$vmat <- vmat.init
param$theta <- vmat.init / rowSums(vmat.init)
param$phi <- phi.init
param$beta <- beta.init
param$lpmu <- logit_pobs_mu.init
param$lpsigma <- logit_pobs_sigma.init
param$pobs <- pobs.init

vmat.post <- theta.post <- theta.jump.post <- theta.accept <- array(, dim=c(nl, nc, niter))
vmat.post[,,1] <- param$vmat
theta.post[,,1] <- param$theta
theta.jump.post[,,1] <- theta.jump
theta.accept[,,1] <- FALSE

phi.post <- phi.jump.post <- phi.accept <- array(, dim=c(nc, ns, niter))
phi.post[,,1] <- param$phi
phi.jump.post[,,1] <- phi.jump
phi.accept[,,1] <- FALSE

lpmu.post <- lpsigma.post <- numeric(niter)
lpmu.post[1] <- param$lpmu
lpsigma.post[1] <- param$lpsigma

pobs.post <- matrix(, ns, niter)
pobs.post[,1] <- param$pobs
pobs.jump.post <- numeric(niter)
pobs.jump.post[1] <- pobs.jump
pobs.accept <- logical(niter)

beta.post <- array(, dim=c(nx,nc,niter))
beta.post[,,1] <- param$beta
beta.jump.post <- numeric(niter)
beta.jump.post[1] <- beta.jump
beta.accept <- logical(niter)

for (i in 2:niter) {
  theta.up <- update.theta(param, jump=theta.jump, nl, nc, ns, y, x, nmat, id.perf)
  vmat.post[,,i] <- param$vmat <- theta.up$vmat
  theta.post[,,i] <- param$theta <- theta.up$theta
  theta.accept[,,i] <- theta.up$accept

  phi.up <- update.phi(param, jump=phi.jump, nl, nc, ns, y, nmat, id.perf, a.phi=1, b.phi=1)
  phi.post[,,i] <- param$phi <- phi.up$phi
  phi.accept[,,i] <- phi.up$accept

  theta.jump.up <- thetaphi.jumpTune(accept=theta.accept, jump=theta.jump, ni=i, adapt=adapt, low=.3, high=.8)
  theta.jump.post[,,i] <- theta.jump <- theta.jump.up$jump

  phi.jump.up <- thetaphi.jumpTune(accept=phi.accept, jump=phi.jump, ni=i, adapt=adapt, low=.3, high=.8)
  phi.jump.post[,,i] <- phi.jump <- phi.jump.up$jump

  lpmu.up <- update.lpmu(param, jump=.1, ns)
  lpmu.post[i] <- param$lpmu <- lpmu.up$lpmu
  lpsigma.post[i] <- param$lpsigma <- lpmu.up$lpsigma

  pobs.up <- update.pobs(param, jump=pobs.jump, nl, ns, id.perf)
  pobs.post[,i] <- param$pobs <- pobs.up$pobs
  pobs.accept[i] <- pobs.up$accept

  pobs.jump.up <- pobsbeta.jumpTune(accept=pobs.accept, jump=pobs.jump, ni=i, adapt=nadapt, low=.3, high=.8)
  pobs.jump.post[i] <- pobs.jump <- pobs.jump.up$jump

  beta.up <- update.beta(param, jump=beta.jump, nl, nc, x)
  beta.post[,,i] <- param$beta <- beta.up$beta
  beta.accept[i] <- beta.up$accept

  beta.jump.up <- pobsbeta.jumpTune(accept=beta.accept, jump=beta.jump, ni=i, adapt=nadapt, low=.3, high=.8)
  beta.jump.post[i] <- beta.jump <- beta.jump.up$jump
}

#==============
# save results
#==============
theta.est <- apply(theta.post[,,nburn:niter], 1:2, median)
phi.est <- apply(phi.post[,,nburn:niter], 1:2, median)
beta.est <- apply(beta.post[,,nburn:niter], 1:2, median)
lpmu.est <- median(lpmu.post[nburn:niter])
lpsigma.est <- median(lpsigma.post[nburn:niter])
pobs.est <- apply(pobs.post[,nburn:niter], 1, median)

res <- list()
res$theta <- theta.est
res$phi <- phi.est
res$beta <- beta.est
res$lpmu <- lpmu.est
res$lpsigma <- lpsigma.est
res$pobs <- pobs.est

class(res) <- c("rlda", "list")
res
} # rlda.regression


