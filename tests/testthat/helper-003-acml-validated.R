## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
## IMPORTANT NOTE
## THIS IS VALIDATED CODE WITH PUBLISHED RESULTS
## DO NOT EDIT, THIS IS A REFERENCE COPY
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ##############################################################################
 #
# ODSMethods Statistical methods in outcome dependent sampling
#
# Copyright (C) 2017 Jonathan Schildcrout
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>

library(MASS)
library(mvtnorm)

av_vi.calc <- function(zi, sigma0, sigma1, rho, sigmae){
  zi %*% matrix(c(sigma0^2, rho*sigma0*sigma1,rho*sigma0*sigma1,sigma1^2), nrow=2) %*% t(zi) +
    sigmae*sigmae*diag(length(zi[,2]))
}
## Ascertainment correction piece for univariate sampling
av_lci <- function(cutpoints, SampProb, mu_q, sigma_q){
  CDFs <- pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
  sum( SampProb*(CDFs[2:length(CDFs)] - CDFs[1:(length(CDFs)-1)]) )
}
## Ascertainment correction piece for bivariate sampling
av_lci.bivar <- function(cutpoints, SampProb, mu_q, sigma_q){
  (SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2]
  #print("blah2")
}
av_subject.ll.lme <- function(yi, xi, beta, vi){
  resid <- yi - xi %*% beta
  -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}

## Calculate log of the ascertainment correction for the univariate sampling case
av_ascertainment.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  vi      <- av_vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- (wi %*% mu)[,1]
  sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
  log(av_lci(cutpoints, SampProb, mu_q, sigma_q))
}
## Calculate log of the ascertainment correction for the bivariate sampling case
av_ascertainment.correction.bivar <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  vi      <- av_vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- as.vector(wi %*% mu)
  sigma_q <- wi %*% vi %*% t(wi)
  log((SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]),
                                        mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2])
}

## Calculate conditional likelihood for the univariate and bivariate sampling cases
av_total.nll.lme <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi){
  total <- 0
  for(i in unique(id)){
    yi <- y[id==i]
    ni <- length(yi)
    xi <- x[id==i,]
    zi <- z[id==i,]
    if (w.function != "bivar"){
      if (w.function=="mean")      wi <- t(rep(1/ni, ni))
      if (w.function=="intercept") wi<- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[1,]
      if (w.function=="slope")     wi<- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[2,]
      wi    <- matrix(wi, 1, ni)
      IPWi  <- 1/ unique(SampProbi[id==i])
      vi    <- av_vi.calc(zi, sigma0, sigma1, rho, sigmae)
      total <- total + av_subject.ll.lme(yi, xi, beta, vi)*IPWi -
        av_ascertainment.correction(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
    }else{
      wi    <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))
      IPWi  <- 1/ unique(SampProbi[id==i])
      vi    <- av_vi.calc(zi, sigma0, sigma1, rho, sigmae)
      total <- total + av_subject.ll.lme(yi, xi, beta, vi)*IPWi -
        av_ascertainment.correction.bivar(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
    }
  }
  -total
}

av_ascertainment.gradient.correction.bivar <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  eps     <- 1e-6
  param   <- c(beta, sigma0, sigma1, rho, sigmae)
  npar    <- length(param)
  vi      <- av_vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- as.vector(wi %*% mu)
  sigma_q <- wi %*% vi %*% t(wi)
  start   <- pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]]

  ## for this bivariate sampling case, right now we calculate gradients numerically
  new.area <- NULL
  eps.mtx <- diag(c(rep(eps,npar)))
  for (rr in 1:npar){
    par.new      <- param+eps.mtx[rr,]
    vi.tmp       <- av_vi.calc(zi, par.new[(npar-3)], par.new[(npar-2)], par.new[(npar-1)], par.new[npar])
    mu.tmp       <- xi %*% par.new[1:(npar-4)]
    mu_q.tmp     <- as.vector(wi %*% mu.tmp)
    sigma_q.tmp  <- wi %*% vi.tmp %*% t(wi)
    new.area <- c(new.area, pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q.tmp, sigma=sigma_q.tmp)[[1]])
  }
  Deriv <- (new.area-start)/eps
  out <- (SampProb[1]-SampProb[2])*Deriv / av_lci.bivar(cutpoints, SampProb, mu_q, sigma_q)
  out
}

av_ascertainment.gradient.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  vi      <- av_vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- (wi %*% mu)[,1]
  sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])

  l <- av_lci(cutpoints, SampProb, mu_q, sigma_q)
  p <- SampProb[1:(length(SampProb)-1)] - SampProb[2:(length(SampProb))]
  f <- dnorm(cutpoints, mu_q, sigma_q)

  d_li_beta <- (wi %*% xi) * sum(p*f) / l

  f_alpha_k <- sum(p*f*(cutpoints - mu_q)) / (l * sigma_q *2 * sqrt(wi %*% vi %*% t(wi))[1,1])
  a5        <- (wi %*% zi %*% matrix(c(2*sigma0, rho*sigma1,    rho*sigma1,     0),        nrow=2) %*% t(zi) %*% t(wi))[1,1]
  a6        <- (wi %*% zi %*% matrix(c(0,        rho*sigma0,    rho*sigma0,     2*sigma1), nrow=2) %*% t(zi) %*% t(wi))[1,1]
  a7        <- (wi %*% zi %*% matrix(c(0,        sigma0*sigma1, sigma0*sigma1,  0),        nrow=2) %*% t(zi) %*% t(wi))[1,1]
  a8        <- (wi %*% (2 * sigmae * diag(length(yi))) %*% t(wi))[1,1]
  c(d_li_beta, c(f_alpha_k * c(a5, a6, a7, a8)))
}

av_subject.gradient.ll.lme <- function(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae){
  resid <- yi - xi %*% beta
  vi    <- av_vi.calc(zi, sigma0, sigma1, rho, sigmae)
  inv.v <- solve(vi)

  dvk5 <- zi %*% matrix(c(2*sigma0, rho*sigma1, rho*sigma1, 0),nrow=2) %*% t(zi)
  dvk6 <- zi %*% matrix(c(0, rho*sigma0, rho*sigma0, 2*sigma1),nrow=2) %*% t(zi)
  dvk7 <- zi %*% matrix(c(0, sigma0*sigma1, sigma0*sigma1, 0),nrow=2) %*% t(zi)
  dvk8 <- 2 * sigmae  * diag(length(yi))

  l14 <- t(xi) %*% inv.v %*% resid # for Beta
  l5  <- -0.5*(sum(diag(inv.v %*% dvk5)) - t(resid) %*% inv.v %*% dvk5 %*% inv.v %*% resid)[1,1]
  l6  <- -0.5*(sum(diag(inv.v %*% dvk6)) - t(resid) %*% inv.v %*% dvk6 %*% inv.v %*% resid)[1,1]
  l7  <- -0.5*(sum(diag(inv.v %*% dvk7)) - t(resid) %*% inv.v %*% dvk7 %*% inv.v %*% resid)[1,1]
  l8  <- -0.5*(sum(diag(inv.v %*% dvk8)) - t(resid) %*% inv.v %*% dvk8 %*% inv.v %*% resid)[1,1]
  list(gr=append(l14, c(l5,l6,l7,l8)),
       vi=vi)
}

#

av_gradient.nll.lme <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=FALSE){
  total <- 0
  param.vec <- c(beta, log(sigma0),log(sigma1),log((1+rho)/(1-rho)),log(sigmae))
  n.par <- length(param.vec)
  cheese <- matrix(0,n.par,n.par)
  for(i in unique(id)){
    yi <- y[id==i]
    ni <- length(yi)
    xi <- x[id==i,]
    zi <- z[id==i,]
    IPWi <- 1/ unique(SampProbi[id==i])
    if (w.function != "bivar"){
      if (w.function=="mean")      wi <- t(rep(1/ni, ni))
      if (w.function=="intercept") wi<- (solve(t(xi[,1:2])%*% xi[,1:2]) %*% t(xi[,1:2]))[1,]
      if (w.function=="slope")     wi<- (solve(t(xi[,1:2])%*% xi[,1:2]) %*% t(xi[,1:2]))[2,]
      wi   <- matrix(wi, 1, ni)
      subject <- av_subject.gradient.ll.lme(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae)
      correct <- av_ascertainment.gradient.correction(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
      Gradi  <- subject[['gr']]*IPWi  + correct ## Gradient for ith subject
      total  <- total + Gradi
    }else{
      wi   <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))
      subject <- av_subject.gradient.ll.lme(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae)
      correct <- av_ascertainment.gradient.correction.bivar(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
      Gradi  <- subject[['gr']]*IPWi  - correct ## Gradient for ith subject: Notice the minus here versus the plus in the univariate case
      total  <- total + Gradi
    }

    #       if(CheeseCalc==TRUE)
    #       {
    #           print(paste("subject",i, "total",total[1]))
    #       }

    if (CheeseCalc==TRUE){
      ## Need to use the chain rule: note that param,vec is on the unconstrained scale but Gradi was calculated on the constrained parameters
      Gradi[(n.par-3)] <- Gradi[(n.par-3)]*exp(param.vec[(n.par-3)])
      Gradi[(n.par-2)] <- Gradi[(n.par-2)]*exp(param.vec[(n.par-2)])
      Gradi[(n.par-1)] <- Gradi[(n.par-1)]*2*exp(param.vec[(n.par-1)])/((exp(param.vec[(n.par-1)])+1)^2)
      Gradi[n.par]     <- Gradi[n.par]  *exp(param.vec[n.par])
      cheese  <- cheese + outer(Gradi, Gradi)
    }
  }
  if (CheeseCalc==TRUE) total <- -cheese
  -total
}


av_LogLikeAndScore <- function(params, y, x, z, id, w.function, cutpoints, SampProb, SampProbi, ProfileCol=NA){
  npar   <- length(params)
  beta   <- params[1:(npar-4)]
  sigma0 <- exp(params[(npar-3)])
  sigma1 <- exp(params[(npar-2)])
  rho    <- (exp(params[(npar-1)])-1) / (exp(params[(npar-1)])+1)
  sigmae <- exp(params[npar])

  out     <- av_total.nll.lme(   y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi)
  GRAD    <- av_gradient.nll.lme(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi)

  ## Need to use the chain rule: note that params is on the unconstrained scale but GRAD was calculated on the constrained parameters
  GRAD[(npar-3)] <- GRAD[(npar-3)]*exp(params[(npar-3)])
  GRAD[(npar-2)] <- GRAD[(npar-2)]*exp(params[(npar-2)])
  GRAD[(npar-1)] <- GRAD[(npar-1)]*2*exp(params[(npar-1)])/((exp(params[(npar-1)])+1)^2)
  GRAD[npar]     <- GRAD[npar]*exp(params[npar])
  ## Force the gradient of the fixed parameter to be zero, so that it does not move
  #if (!is.na(ProfileCol)) GRAD[ProfileCol] <- 0
  attr(out,"gradient") <- GRAD

  out
}

## If you do not want to use the ascertainment correction term in the conditional likelihood
## set all SampProb values equal to each other.  This would be the case if you were doing
## straightforward maximum likelihood (albeit inefficient) or weighted likelihood.
acml_validated <- function(
  y,                           ## response vector
  x,                           ## fixed effects design matrix
  z,                           ## random effects design matrix (right now this should be an intercept and a time-varying covariate (?)
  id,                          ## subject id variable
  w.function="mean",           ## Function upon which sampling is based. Choices are the univariate "mean", "intercept", "slope", and the bivariate "bivar"
  InitVals,                    ## Starting values
  cutpoints = c(0,5),          ## the cutpoints: when w.function="bivar", this is a vector of length 4 that define a central, rectangular region with vertices (x_lower, x_upper, y_lower, y_upper).
  SampProb = c(1, 1, 1),       ## Sampling probabilities within the sampling strata to be used for ACML
  SampProbi=rep(1, length(y)), ## Subject specific sampling probabilities to only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar
  ProfileCol=NA)               ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
{
  out <- nlm(av_LogLikeAndScore, InitVals, y=y, x=x, z=z, id=id, w.function=w.function, cutpoints=cutpoints, SampProb=SampProb,
             SampProbi=SampProbi, ProfileCol=ProfileCol, gradtol=1e-12,
             stepmax=4, iterlim=250, check.analyticals = TRUE) #, print.level=1)

  ## Calculate the observed information and then invert to get the covariance matrix
  npar <- length(out$estimate)
  Hessian.eps <- 1e-7
  eps.mtx     <- diag(rep(Hessian.eps, npar))
  grad.at.max <- out$gradient
  ObsInfo     <- matrix(NA, npar, npar)

  ## Observed Information
  for (j in 1:npar){
    temp        <- av_LogLikeAndScore(out$estimate+eps.mtx[j,], y=y, x=x, z=z, id=id,
                                   w.function=w.function, cutpoints=cutpoints,
                                   SampProb=SampProb,SampProbi=SampProbi, ProfileCol=ProfileCol)
    ObsInfo[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
  }
  ## Cheese part of the sandwich estimator
  Cheese <- av_gradient.nll.lme(y=y, x=x, z=z, w.function=w.function,
                             id=id, beta=out$estimate[c(1:(npar-4))],
                             sigma0=exp(out$estimate[(npar-3)]),
                             sigma1=exp(out$estimate[(npar-2)]),
                             rho=   (exp(out$estimate[(npar-1)])-1) / (exp(out$estimate[(npar-1)])+1),
                             sigmae=exp(out$estimate[npar]),
                             cutpoints=cutpoints,
                             SampProb=SampProb,
                             SampProbi=SampProbi,
                             CheeseCalc=TRUE)

  if (!is.na(ProfileCol)){
    out$estimate <- out$estimate[-ProfileCol]
    ObsInfo <- ObsInfo[-ProfileCol, -ProfileCol]
    Cheese  <- Cheese[-ProfileCol, -ProfileCol]
  }

  list(Ests=out$estimate, covar=solve(ObsInfo), LogL= -out$minimum,  Code=out$code, robcov=solve(ObsInfo)%*%Cheese%*%solve(ObsInfo))
}

## BDS specific validation functions. These reuse the common ACML validation functions above and add the phase-I BLUP quantities used by BDS.

av_bds_logACi1q <- function(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e,
                          cutpoints, SampProb, mui.phase1)
{
  vi <- av_vi.calc(zi, sigma.vc[1], sigma.vc[2], rho.vc[1], sigma.e[1])
  mu <- xi %*% beta
  mu_q <- (wi %*% (mu - mui.phase1))[, 1]
  sigma_q <- sqrt((wi %*% vi %*% t(wi))[1, 1])

  list(
    vi = vi,
    logACi = log(av_lci(cutpoints, SampProb, mu_q, sigma_q))
  )
}

av_bds_LogLikeiC2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e)
{
  yi <- subjectData[["yi"]]
  xi <- subjectData[["xi"]]
  zi <- subjectData[["zi"]]
  wi <- subjectData[["wi"]]
  SampProb <- subjectData[["SampProb.i"]]
  cutpoints <- subjectData[["cutpoints.i"]]
  mui.phase1 <- subjectData[["mui.phase1"]]
  Weights.i <- subjectData[["Weights.i"]]

  tmp <- av_bds_logACi1q(
    yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e,
    cutpoints, SampProb, mui.phase1
  )
  logACi <- tmp[["logACi"]]
  liC <- av_subject.ll.lme(yi, xi, beta, tmp[["vi"]]) * Weights.i - logACi
  c(liC, logACi)
}

av_bds_LogLikeC2 <- function(y, x, z, w.function, id, beta, sigma.vc, rho.vc,
                           sigma.e, cutpoints, SampProb, Weights,
                           Keep.liC = FALSE, xcol.phase1, ests.phase1)
{
  subjectData <- av_bds_CreateSubjectData(
    id = id, y = y, x = x, z = z, Weights = Weights, SampProb = SampProb,
    cutpoints = cutpoints, w.function = w.function,
    xcol.phase1 = xcol.phase1, ests.phase1 = ests.phase1
  )
  liC.and.logACi <- lapply(
    subjectData, av_bds_LogLikeiC2,
    beta = beta, sigma.vc = sigma.vc, rho.vc = rho.vc, sigma.e = sigma.e
  )

  if (Keep.liC == FALSE) {
    -1 * Reduce("+", liC.and.logACi)[1]
  } else {
    list(
      liC = c(unlist(sapply(liC.and.logACi, function(x) x[1]))),
      logACi = c(unlist(sapply(liC.and.logACi, function(x) x[2])))
    )
  }
}

av_bds_logACi1q.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e)
{
  xi <- subjectData[["xi"]]
  zi <- subjectData[["zi"]]
  wi <- subjectData[["wi"]]
  SampProb <- subjectData[["SampProb.i"]]
  cutpoints <- subjectData[["cutpoints.i"]]
  mui.phase1 <- subjectData[["mui.phase1"]]
  ni <- nrow(zi)

  wi <- matrix(wi, 1, ni)
  vi <- av_vi.calc(zi, sigma.vc[1], sigma.vc[2], rho.vc[1], sigma.e[1])
  t.wi <- t(wi)
  wi.zi <- wi %*% zi
  t.wi.zi <- t(wi.zi)

  mu <- xi %*% beta
  mu_q <- (wi %*% (mu - mui.phase1))[, 1]
  sigma_q <- sqrt((wi %*% vi %*% t.wi)[1, 1])

  l <- av_lci(cutpoints, SampProb, mu_q, sigma_q)
  p <- SampProb[1:(length(SampProb) - 1)] - SampProb[2:(length(SampProb))]
  f <- dnorm(cutpoints, mu_q, sigma_q)

  d_li_beta <- (wi %*% xi) * sum(p * f) / l
  f_alpha_k <- sum(p * f * (cutpoints - mu_q)) / (l * 2 * sigma_q^2)

  len.sigma.vc <- length(sigma.vc)
  len.rho.vc <- length(rho.vc)
  len.sigma.e <- length(sigma.e)

  SDMat.RE <- diag(sigma.vc)
  b0 <- matrix(0, len.sigma.vc, len.sigma.vc)
  b0[lower.tri(b0, diag = FALSE)] <- rho.vc
  CorMat.RE <- t(b0) + b0 + diag(rep(1, len.sigma.vc))
  D <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
  m <- matrix(0, len.sigma.vc, len.sigma.vc)

  dVi.dsigma.vc <- NULL
  for (mmm in 1:len.sigma.vc) {
    m1 <- m
    m1[mmm, ] <- m1[, mmm] <- 1
    m1[mmm, mmm] <- 2
    tmp <- sigma.vc[mmm] + 1 * (sigma.vc[mmm] == 0)
    dViMat.dsigma.vc <- m1 * D / tmp
    dVi.dsigma.vc <- c(
      dVi.dsigma.vc,
      (wi.zi %*% dViMat.dsigma.vc %*% t.wi.zi)[1, 1]
    )
  }

  b1 <- matrix(0, len.sigma.vc, len.sigma.vc)
  b1[lower.tri(b1, diag = FALSE)] <- c(1:len.rho.vc)
  tbb1 <- t(b1) + b1

  dVi.drho.vc <- NULL
  for (mmm in 1:len.rho.vc) {
    tmp <- which(tbb1 == mmm, arr.ind = TRUE)
    m2 <- m
    m2[tmp[1, 1], tmp[1, 2]] <- 1
    m2[tmp[2, 1], tmp[2, 2]] <- 1
    tmp <- rho.vc[mmm] + 1 * (rho.vc[mmm] == 0)
    dViMat.drho.vc <- m2 * D / tmp
    dVi.drho.vc <- c(
      dVi.drho.vc,
      (wi.zi %*% dViMat.drho.vc %*% t.wi.zi)[1, 1]
    )
  }

  dVi.dsigma.e <- NULL
  for (mmm in 1:len.sigma.e) {
    dsigma.e.vec <- 2 * sigma.e
    dsigma.e.vec[-mmm] <- 0
    dViMat.dsigma.e <- diag(rep(dsigma.e.vec, each = ni / len.sigma.e))
    dVi.dsigma.e <- c(
      dVi.dsigma.e,
      (wi %*% dViMat.dsigma.e %*% t.wi)[1, 1]
    )
  }

  c(
    d_li_beta,
    f_alpha_k * c(dVi.dsigma.vc, dVi.drho.vc, dVi.dsigma.e)
  )
}

av_bds_li.lme.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e)
{
  score <- av_subject.gradient.ll.lme(
    yi = subjectData[["yi"]],
    xi = subjectData[["xi"]],
    zi = subjectData[["zi"]],
    beta = beta,
    sigma0 = sigma.vc[1],
    sigma1 = sigma.vc[2],
    rho = rho.vc[1],
    sigmae = sigma.e[1]
  )
  score[["gr"]] <- score[["gr"]] * subjectData[["Weights.i"]]
  score
}

av_bds_LogLikeC.Score2 <- function(y, x, z, w.function, id, beta, sigma.vc,
                                 rho.vc, sigma.e, cutpoints, SampProb,
                                 Weights, CheeseCalc = FALSE,
                                 xcol.phase1, ests.phase1)
{
  param.vec <- c(beta, log(sigma.vc), log((1 + rho.vc) / (1 - rho.vc)), log(sigma.e))
  npar <- length(param.vec)

  len.beta <- length(beta)
  len.sigma.vc <- length(sigma.vc)
  len.rho.vc <- length(rho.vc)
  len.sigma.e <- length(sigma.e)

  vc.sd.index <- len.beta + (c(1:len.sigma.vc))
  vc.rho.index <- len.beta + len.sigma.vc + (c(1:len.rho.vc))
  err.sd.index <- len.beta + len.sigma.vc + len.rho.vc + c(1:len.sigma.e)
  notbeta.index <- c(vc.sd.index, vc.rho.index, err.sd.index)

  subjectData <- av_bds_CreateSubjectData(
    id = id, y = y, x = x, z = z, Weights = Weights, SampProb = SampProb,
    cutpoints = cutpoints, w.function = w.function,
    xcol.phase1 = xcol.phase1, ests.phase1 = ests.phase1
  )

  UncorrectedScorei <- lapply(
    subjectData, av_bds_li.lme.score2,
    beta = beta, sigma.vc = sigma.vc, rho.vc = rho.vc, sigma.e = sigma.e
  )
  Gradienti <- lapply(UncorrectedScorei, function(x) x[["gr"]])
  UncorrectedScore <- Reduce("+", Gradienti)

  logACi.Score <- lapply(
    subjectData, av_bds_logACi1q.score2,
    beta = beta, sigma.vc = sigma.vc, rho.vc = rho.vc, sigma.e = sigma.e
  )
  logAC.Score <- Reduce("+", logACi.Score)
  CorrectedScore <- UncorrectedScore + logAC.Score

  if (CheeseCalc == TRUE) {
    GradiMat <- mapply("+", Gradienti, logACi.Score)
    GradiMat[notbeta.index, ] <- GradiMat[notbeta.index, ] *
      c(exp(param.vec[vc.sd.index]),
        2 * exp(param.vec[vc.rho.index]) / ((exp(param.vec[vc.rho.index]) + 1)^2),
        exp(param.vec[err.sd.index]))
    cheese <- matrix(0, npar, npar)
    for (mm in 1:ncol(GradiMat)) {
      cheese <- cheese + outer(GradiMat[, mm], GradiMat[, mm])
    }
  }

  out <- -CorrectedScore
  if (CheeseCalc == TRUE) out <- cheese
  out
}

av_bds_LogLikeCAndScore2 <- function(params, y, x, z, id, w.function, cutpoints,
                                   SampProb, Weights, ProfileCol = NA,
                                   Keep.liC = FALSE, xcol.phase1, ests.phase1)
{
  npar <- length(params)

  nbeta <- ncol(x)
  nVCsd <- ncol(z)
  nVCrho <- choose(nVCsd, 2)
  nERRsd <- npar - nbeta - nVCsd - nVCrho

  beta.index <- c(1:nbeta)
  vc.sd.index <- nbeta + (c(1:nVCsd))
  vc.rho.index <- nbeta + nVCsd + (c(1:nVCrho))
  err.sd.index <- nbeta + nVCsd + nVCrho + c(1:nERRsd)

  beta <- params[beta.index]
  sigma.vc <- exp(params[vc.sd.index])
  rho.vc <- (exp(params[vc.rho.index]) - 1) / (exp(params[vc.rho.index]) + 1)
  sigma.e <- exp(params[err.sd.index])

  out <- av_bds_LogLikeC2(
    y = y, x = x, z = z, w.function = w.function, id = id,
    beta = beta, sigma.vc = sigma.vc, rho.vc = rho.vc, sigma.e = sigma.e,
    cutpoints = cutpoints, SampProb = SampProb, Weights = Weights,
    Keep.liC = Keep.liC, xcol.phase1 = xcol.phase1, ests.phase1 = ests.phase1
  )
  GRAD <- av_bds_LogLikeC.Score2(
    y = y, x = x, z = z, w.function = w.function, id = id,
    beta = beta, sigma.vc = sigma.vc, rho.vc = rho.vc, sigma.e = sigma.e,
    cutpoints = cutpoints, SampProb = SampProb, Weights = Weights,
    xcol.phase1 = xcol.phase1, ests.phase1 = ests.phase1
  )

  GRAD[vc.sd.index] <- GRAD[vc.sd.index] * exp(params[vc.sd.index])
  GRAD[vc.rho.index] <- GRAD[vc.rho.index] *
    2 * exp(params[vc.rho.index]) / ((exp(params[vc.rho.index]) + 1)^2)
  GRAD[err.sd.index] <- GRAD[err.sd.index] * exp(params[err.sd.index])
  if (!is.na(ProfileCol)) GRAD[ProfileCol] <- 0
  attr(out, "gradient") <- GRAD

  out
}

av_bds_CreateSubjectData <- function(id, y, x, z, Weights, SampProb, cutpoints,
                                   w.function, xcol.phase1, ests.phase1)
{
  npar.phase1 <- length(ests.phase1)
  x.phase1 <- x[, xcol.phase1]
  nbeta.phase1 <- ncol(x.phase1)
  nVCsd.phase1 <- ncol(z)
  nVCrho.phase1 <- choose(nVCsd.phase1, 2)
  nERRsd.phase1 <- npar.phase1 - nbeta.phase1 - nVCsd.phase1 - nVCrho.phase1

  beta.index.phase1 <- c(1:nbeta.phase1)
  vc.sd.index.phase1 <- nbeta.phase1 + (c(1:nVCsd.phase1))
  vc.rho.index.phase1 <- nbeta.phase1 + nVCsd.phase1 + (c(1:nVCrho.phase1))
  err.sd.index.phase1 <- nbeta.phase1 + nVCsd.phase1 + nVCrho.phase1 + c(1:nERRsd.phase1)

  beta.phase1 <- ests.phase1[beta.index.phase1]
  sigma.vc.phase1 <- exp(ests.phase1[vc.sd.index.phase1])
  rho.vc.phase1 <- (exp(ests.phase1[vc.rho.index.phase1]) - 1) /
    (exp(ests.phase1[vc.rho.index.phase1]) + 1)
  sigma.e.phase1 <- exp(ests.phase1[err.sd.index.phase1])

  D.tmp <- diag(sigma.vc.phase1)
  b <- matrix(0, nVCsd.phase1, nVCsd.phase1)
  b[lower.tri(b, diag = FALSE)] <- rho.vc.phase1
  R.tmp <- t(b) + b + diag(rep(1, nVCsd.phase1))
  D.phase1 <- D.tmp %*% R.tmp %*% D.tmp
  x.phase1.tmp <- split(x.phase1, id)
  ncol.x.phase1 <- ncol(x.phase1)

  id.tmp <- split(id, id)
  y.tmp <- split(y, id)
  x.tmp <- split(x, id)
  z.tmp <- split(z, id)
  Weights.tmp <- split(Weights, id)
  SampProb.tmp <- split(SampProb, id)
  cutpoints.tmp <- split(cutpoints, id)
  w.function.tmp <- split(w.function, id)

  ncol.x <- ncol(x)
  ncol.z <- ncol(z)
  ncol.SampProb <- ncol(SampProb)
  ncol.cutpoints <- ncol(cutpoints)

  subjectData <- list()
  uid <- as.character(unique(id))
  for (j in seq(along = uid)) {
    i <- uid[j]
    zi <- matrix(z.tmp[[i]], ncol = ncol.z)
    xi <- matrix(x.tmp[[i]], ncol = ncol.x)
    yi <- y.tmp[[i]]
    w.functioni <- unique(w.function.tmp[[i]])
    ni <- length(yi)

    xi.phase1 <- matrix(x.phase1.tmp[[i]], ncol = ncol.x.phase1)
    Vi.phase1 <- av_vi.calc(
      zi = zi,
      sigma0 = sigma.vc.phase1[1],
      sigma1 = sigma.vc.phase1[2],
      rho = rho.vc.phase1[1],
      sigmae = sigma.e.phase1[1]
    )
    wi.tmp <- D.phase1 %*% t(zi) %*% solve(Vi.phase1)
    mui.phase1 <- xi.phase1 %*% beta.phase1

    if (w.functioni == "blup.intercept") {
      wi <- wi.tmp[1, ]
    } else if (w.functioni == "blup.slope") {
      wi <- wi.tmp[2, ]
    } else if (w.functioni == "mean") {
      wi <- t(rep(1 / ni, ni))
      mui.phase1 <- rep(0, nrow(xi))
    } else {
      stop("BDS validation helper only supports mean, blup.intercept, and blup.slope")
    }
    wi <- matrix(wi, 1, ni)

    subjectData[[j]] <- list(
      idi = as.character(unique(id.tmp[[i]])),
      xi = xi,
      zi = zi,
      yi = yi,
      wi = wi,
      mui.phase1 = mui.phase1,
      Weights.i = unique(Weights.tmp[[i]]),
      SampProb.i = matrix(SampProb.tmp[[i]], ncol = ncol.SampProb)[1, ],
      w.function.i = as.character(w.functioni),
      cutpoints.i = matrix(cutpoints.tmp[[i]], ncol = ncol.cutpoints)[1, ]
    )
  }
  names(subjectData) <- uid
  subjectData
}

av_bds_acml.lmem2 <- function(formula.fixed, formula.random, data, id,
                            w.function = "mean", InitVals, cutpoints,
                            SampProb, Weights, ProfileCol = NA,
                            xcol.phase1 = NULL, ests.phase1 = NULL)
{
  terms <- unique(c(
    all.vars(formula.fixed),
    all.vars(formula.random),
    as.character(substitute(id)),
    as.character(substitute(Weights))
  ))
  data <- data[, terms]

  if (any(is.na(data))) data <- na.omit(data)

  id0 <- as.character(substitute(id))
  id <- data$id <- data[, id0]

  fixed.f <- model.frame(formula.fixed, data)
  y <- model.response(fixed.f, "numeric")
  x <- model.matrix(formula.fixed, fixed.f)
  z <- model.matrix(formula.random, fixed.f)

  Weights0 <- as.character(substitute(Weights))
  Weights <- data$Weights <- data[, Weights0]

  acml.fit <- nlm(
    av_bds_LogLikeCAndScore2,
    InitVals,
    y = y,
    x = x,
    z = z,
    id = id,
    w.function = w.function,
    cutpoints = cutpoints,
    SampProb = SampProb,
    Weights = Weights,
    ProfileCol = ProfileCol,
    xcol.phase1 = xcol.phase1,
    ests.phase1 = ests.phase1,
    stepmax = 4,
    iterlim = 250,
    check.analyticals = TRUE,
    print.level = 0
  )

  npar <- length(acml.fit$estimate)
  Hessian.eps <- 1e-7
  eps.mtx <- diag(rep(Hessian.eps, npar))
  grad.at.max <- acml.fit$gradient
  ObsInfo.tmp <- ObsInfo <- matrix(NA, npar, npar)

  for (j in 1:npar) {
    temp <- av_bds_LogLikeCAndScore2(
      acml.fit$estimate + eps.mtx[j, ],
      y = y,
      x = x,
      z = z,
      id = id,
      w.function = w.function,
      cutpoints = cutpoints,
      SampProb = SampProb,
      Weights = Weights,
      ProfileCol = ProfileCol,
      xcol.phase1 = xcol.phase1,
      ests.phase1 = ests.phase1
    )
    ObsInfo.tmp[j, ] <- (attr(temp, "gradient") - grad.at.max) / Hessian.eps
  }
  for (m in 1:npar) {
    for (n in 1:npar) {
      ObsInfo[m, n] <- (ObsInfo.tmp[m, n] + ObsInfo.tmp[n, m]) / 2
    }
  }

  nbeta <- ncol(x)
  nVCsd <- ncol(z)
  nVCrho <- choose(nVCsd, 2)
  nERRsd <- npar - nbeta - nVCsd - nVCrho

  beta.index <- c(1:nbeta)
  vc.sd.index <- nbeta + (c(1:nVCsd))
  vc.rho.index <- nbeta + nVCsd + (c(1:nVCrho))
  err.sd.index <- nbeta + nVCsd + nVCrho + c(1:nERRsd)

  Cheese <- av_bds_LogLikeC.Score2(
    y = y,
    x = x,
    z = z,
    w.function = w.function,
    id = id,
    beta = acml.fit$estimate[beta.index],
    sigma.vc = exp(acml.fit$estimate[vc.sd.index]),
    rho.vc = (exp(acml.fit$estimate[vc.rho.index]) - 1) /
      (exp(acml.fit$estimate[vc.rho.index]) + 1),
    sigma.e = exp(acml.fit$estimate[err.sd.index]),
    cutpoints = cutpoints,
    SampProb = SampProb,
    Weights = Weights,
    CheeseCalc = TRUE,
    xcol.phase1 = xcol.phase1,
    ests.phase1 = ests.phase1
  )

  if (!is.na(ProfileCol)) {
    acml.fit$estimate <- acml.fit$estimate[-ProfileCol]
    ObsInfo <- ObsInfo[-ProfileCol, -ProfileCol]
    Cheese <- Cheese[-ProfileCol, -ProfileCol]
  }

  ObsInfo.inv <- solve(ObsInfo)
  list(
    coefficients = acml.fit$estimate,
    covariance = ObsInfo.inv,
    robcov = ObsInfo.inv %*% Cheese %*% ObsInfo.inv,
    logLik = -acml.fit$minimum,
    Code = acml.fit$code
  )
}
