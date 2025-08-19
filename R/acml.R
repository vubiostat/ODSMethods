  ##############################################################################
 #
# ODSMethods Statistical methods in outcome dependent sampling
#
# Copyright (c) 2017 Jonathan Schildcrout
# Copyright (c) 2025 Shawn P. Garbett, Bailu Yan
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

# NOTE: This is the editable version of ACML fitting. The tests will make
# sure that it returns the right results.

vi.calc <- function(zi, sigma0, sigma1, rho, sigmae){
  zi %*% matrix(c(sigma0^2, rho*sigma0*sigma1,rho*sigma0*sigma1,sigma1^2), nrow=2) %*% t(zi) +
    sigmae*sigmae*diag(length(zi[,2]))
}
## Ascertainment correction piece for univariate sampling
lci <- function(cutpoints, SampProb, mu_q, sigma_q){
  CDFs <- pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
  sum( SampProb*(CDFs[2:length(CDFs)] - CDFs[1:(length(CDFs)-1)]) )
}

## Ascertainment correction piece for bivariate sampling
#' @importFrom mvtnorm pmvnorm
lci.bivar <- function(cutpoints, SampProb, mu_q, sigma_q){
  (SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2]
  #print("blah2")
}
subject.ll.lme <- function(yi, xi, beta, vi){
  resid <- yi - xi %*% beta
  -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}

## Calculate log of the ascertainment correction for the univariate sampling case
ascertainment.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- (wi %*% mu)[,1]
  sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
  log(lci(cutpoints, SampProb, mu_q, sigma_q))
}
## Calculate log of the ascertainment correction for the bivariate sampling case
ascertainment.correction.bivar <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- as.vector(wi %*% mu)
  sigma_q <- wi %*% vi %*% t(wi)
  log((SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]),
                                        mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2])
}

## Calculate conditional likelihood for the univariate and bivariate sampling cases
total.nll.lme <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi){
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
      vi    <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
      total <- total + subject.ll.lme(yi, xi, beta, vi)*IPWi -
        ascertainment.correction(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
    }else{
      wi    <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))
      IPWi  <- 1/ unique(SampProbi[id==i])
      vi    <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
      total <- total + subject.ll.lme(yi, xi, beta, vi)*IPWi -
        ascertainment.correction.bivar(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
    }
  }
  -total
}

ascertainment.gradient.correction.bivar <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  eps     <- 1e-6
  param   <- c(beta, sigma0, sigma1, rho, sigmae)
  npar    <- length(param)
  vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- as.vector(wi %*% mu)
  sigma_q <- wi %*% vi %*% t(wi)
  start   <- pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]]

  ## for this bivariate sampling case, right now we calculate gradients numerically
  new.area <- NULL
  eps.mtx <- diag(c(rep(eps,npar)))
  for (rr in 1:npar){
    par.new      <- param+eps.mtx[rr,]
    vi.tmp       <- vi.calc(zi, par.new[(npar-3)], par.new[(npar-2)], par.new[(npar-1)], par.new[npar])
    mu.tmp       <- xi %*% par.new[1:(npar-4)]
    mu_q.tmp     <- as.vector(wi %*% mu.tmp)
    sigma_q.tmp  <- wi %*% vi.tmp %*% t(wi)
    new.area <- c(new.area, pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q.tmp, sigma=sigma_q.tmp)[[1]])
  }
  Deriv <- (new.area-start)/eps
  out <- (SampProb[1]-SampProb[2])*Deriv / lci.bivar(cutpoints, SampProb, mu_q, sigma_q)
  out
}

#' @importFrom stats dnorm pnorm
ascertainment.gradient.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
  vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
  mu      <- xi %*% beta
  mu_q    <- (wi %*% mu)[,1]
  sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])

  l <- lci(cutpoints, SampProb, mu_q, sigma_q)
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

subject.gradient.ll.lme <- function(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae){
  resid <- yi - xi %*% beta
  vi    <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
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

gradient.nll.lme <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=FALSE){
  total <- 0
  param.vec <- c(beta, log(sigma0),log(sigma1),2*atanh(rho),log(sigmae))
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
      subject <- subject.gradient.ll.lme(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae)
      correct <- ascertainment.gradient.correction(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
      Gradi  <- subject[['gr']]*IPWi  + correct ## Gradient for ith subject
      total  <- total + Gradi
    }else{
      wi   <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))
      subject <- subject.gradient.ll.lme(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae)
      correct <- ascertainment.gradient.correction.bivar(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
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


LogLikeAndScore <- function(params, y, x, z, id, w.function, cutpoints, SampProb, SampProbi, ProfileCol=NA){
  npar   <- length(params)
  beta   <- params[1:(npar-4)]
  sigma0 <- exp(params[(npar-3)])
  sigma1 <- exp(params[(npar-2)])
  rho    <- (exp(params[(npar-1)])-1) / (exp(params[(npar-1)])+1)
  sigmae <- exp(params[npar])

  out     <- total.nll.lme(   y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi)
  GRAD    <- gradient.nll.lme(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi)

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
#' @importFrom stats nlm
acml_internal <- function(
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
  out <- nlm(LogLikeAndScore, InitVals, y=y, x=x, z=z, id=id, w.function=w.function, cutpoints=cutpoints, SampProb=SampProb,
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
    temp        <- LogLikeAndScore(out$estimate+eps.mtx[j,], y=y, x=x, z=z, id=id,
                                   w.function=w.function, cutpoints=cutpoints,
                                   SampProb=SampProb,SampProbi=SampProbi, ProfileCol=ProfileCol)
    ObsInfo[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
  }
  ## Cheese part of the sandwich estimator
  Cheese <- gradient.nll.lme(y=y, x=x, z=z, w.function=w.function,
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
  Ests = out$estimate# Lucy updated on 5/30/2025
  names(Ests) = c(colnames(x), "log(sigma_0)","log(sigma_1)","2*atanh(rho)","log(sigmae)")
  list(Ests=Ests, covar=solve(ObsInfo), LogL= -out$minimum,  Code=out$code, robcov=solve(ObsInfo)%*%Cheese%*%solve(ObsInfo))
}

  ##############################################################################
 #
# S3 methods for acml

triangle <- function(n) n*(n+1)/2

#' @export
#' @rdname coef
fixef <- function(object, ...) UseMethod("fixef")

#' @export
fixef.acml <- function(object, ...)
{
  est <- object$Est
  est[1:(length(est)-triangle(object$design$n_rand)-1)]
}

ranef_transform <- function(ran, n_rand)
{
  x        <- exp(ran)
  rng      <- (n_rand+1):(length(ran) - 1)
  x[rng]   <- tanh(ran[rng]/2)
  names(x) <- gsub(".*\\((.*)\\).*", "\\1", names(x), perl=TRUE)
  x
}

#' @export
#' @rdname coef
ranef <- function(object, transform=FALSE, ...) UseMethod("ranef")

#' @export
ranef.acml <- function(object, transform=FALSE, ...)
{
  est <- object$Est
  le  <- length(est)
  ran <- est[(le - triangle(object$design$n_rand)):le]

  if(transform) ranef_transform(ran, object$design$n_rand) else ran
}

#' Extract parameters from fitted models
#'
#' These are S3 methods to extract the entire parameter set, just the fixed
#' effects, or just the random effects. They are by default returned on the
#' unconstrained optimization scale.
#'
#' @exportS3Method
#' @rdname coef
#' @param object the fitted model object to extract model coefficients.
#' @param complete Not used at present, required for S3
#' @param transform logical(1); If TRUE the coefficients will be inverse
#' transformed back to their original scale, this is `exp` for deviation
#' components and Fisher transformed for correlation.
#' @param ... Additional arguments passed along
#' @return A named vector of desired coeffients.
#' @examples
#' data(gbti)
#' design <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
#'               data=gbti, quantiles=c(0.1, 0.9))
#' est <- acml(Response ~ Month*Genotype, design, gbti)
#' coef(est)
#' ranef(est)
#' fixef(est)
coef.acml <- function(object, complete = TRUE, transform = FALSE, ...)
{
  c(fixef(object, ...),
    ranef(object, transform=transform, ...))
}

#' @exportS3Method
vcov.acml <- function(object, complete = TRUE, robust = FALSE, ...)
{
  nm <- names(coef(object))
  vc <- if(robust) object$robcov else object$covar
  rownames(vc) <- nm
  colnames(vc) <- nm
  vc
}

#' @exportS3Method
print.acml <- function(x, digits = max(3L, getOption("digits")), transform = FALSE, ...)
{
  object <- x
  cat("\nCalls:\n",
      paste(deparse(object$design$call), collapse="\n"),
      "\n",
      paste(deparse(object$call), collapse="\n"),
      "\n\n",
      "Cutpoints:\n",
      sep="")
  print(round(object$design$cutpoints, digits=digits), ...)
  cat("\nFixed Effects:\n")
  print(round(fixef(object), digits=digits), ...)
  cat("\nRandom Effects:\n")
  print(round(ranef(object, transform=transform), digits=digits), ...)
  cat("\nNumber of Subjects:\n")
  print(length(unique(object$design$model.frame[,object$design$id])))
  invisible(object)
}

#' @exportS3Method
#' @importFrom stats qnorm vcov
summary.acml <- function(object, digits = max(3L, getOption("digits")),
                         transform = FALSE, robust = FALSE, ...)
{
  raw              <- coef(object, ...)
  beta             <- coef(object, transform=transform, ...)
  se               <- sqrt(diag(vcov(object)))
  object$transform <- transform
  object$robust    <- robust
  object$digits    <- digits
  z                <- raw/se
  le               <- length(beta)
  object$n_random  <- triangle(object$design$n_rand)+1
  object$n_fixed   <- le - object$n_random

  object$coefficients  <- cbind(
    Estimate     = beta,
    `Std. Error` = se,
    L95          = raw + se*qnorm(0.025),
    U95          = raw + se*qnorm(0.975),
    `z value`    = z,
    `Pr(>|z|)`   = 2 * pnorm(abs(z), lower.tail=FALSE)
  )

  if(transform)
  {
    n_rand <- object$design$n_rand

    ran <- (le - triangle(n_rand)):nrow(object$coefficients)

    object$coefficients[ran, 'L95'] <- ranef_transform(object$coefficients[ran, 'L95'], n_rand)
    object$coefficients[ran, 'U95'] <- ranef_transform(object$coefficients[ran, 'U95'], n_rand)

    # Delta Method SE[f(x_hat)] ~= |g'(x_hat)| * SE[x_hat]

    # The one with special handling is the rho(s)
    rhos <- ran[(n_rand+1):(length(ran) - 1)] # Find rho(s)
    object$coefficients[rhos, 2] <- abs(1/cosh(coef(object)[rhos]/2)^2/2) *object$coefficients[rhos, 2]

    # This works for rest coefficients since they are just exp(x) and D(exp(x)) = exp(x)
    ran <- ran[!ran %in% rhos]
    object$coefficients[ran, 2] <- object$coefficients[ran, 1]*object$coefficients[ran, 2]
  }

  # object$residuals <- residuals(object)


  class(object)    <- c("summary.acml", "acml")
  object
}

#' @exportS3Method
#' @importFrom stats printCoefmat
print.summary.acml <- function(x, digits=NULL, signif.stars = getOption("show.signif.stars"), ...)
{
  object <- x
  if(is.null(digits)) digits <- object$digits

  cat("\nCalls:\n",
      paste(deparse(object$design$call), collapse="\n"),
      "\n",
      paste(deparse(object$call), collapse="\n"),
      "\n\n",
      "Cutpoints:\n",
      sep="")
  print(round(object$design$cutpoints, digits=digits), ...)
  cat("\nFixed Effects:\n")
  printCoefmat(object$coefficients[1:object$n_fixed,], digits = digits-1, dig.tst=digits, signif.stars = signif.stars,
               na.print = "NA", ...)

  cat("\nRandom Effects:\n")
  random <- object$coefficients[(object$n_fixed+1):nrow(object$coefficients),]
  nt     <- max(nchar(object$design$id)+1, 7) # Groups and a space is minimum
  pad <- paste0(rep(" ", nt), collapse="")
  rownames(random) <-
    if(object$transform)
    {
      c(
        paste(object$design$id, "(Intercept)"),
        paste0(pad, object$design$time),
        paste0(pad, "rho"),
        "Residual"
      )
    } else
    {
      c(
        paste(object$design$id, "log(Intercept)"),
        paste0(pad, "log(", object$design$time, ")"),
        paste0(pad, "2*atanh(rho)"),
        "log(Residual)"
      )
    }

  print(round(random[,1:4], digits))

  cat("\nNumber of Subjects:", length(unique(object$design$model.frame[,object$design$id])))
  cat("\nNumber of Observations:", nrow(object$design$model.frame))
  cat("\n")
  if(object$transform) cat("\nStd. errors approximated via delta method. Confidence intervals are computed on transformed scale and transformed back and will not be symmetric.\n\n")
  invisible(object)
}

rand.effect.matrix <- function(gamma)
{
  sigma_0 <- gamma[1]
  sigma_1 <- gamma[2]
  rho     <- gamma[3]

  diag(c(sigma_0, sigma_1)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(c(sigma_0, sigma_1))
}

#' @exportS3Method
#' @importFrom stats na.action predict
predict.acml <- function(object, digits=NULL,  ...)
{
  x_mm           <- object$model.matrix
  y              <- object$response
  beta           <- fixef(object)
  y_hat_fixed    <- x_mm %*% matrix(beta, ncol = 1)
  gamma          <- ranef(object, transform=TRUE)
  G              <- rand.effect.matrix(gamma)
  sigmae2        <- gamma[length(gamma)]^2
  subject_ids    <- object$ids
  n_subjects     <- length(unique(subject_ids))
  y_hat_random   <- matrix(NA, nrow = nrow(x_mm), ncol = 1)

  for (j in 1:n_subjects)
  {
    subject_id = unique(subject_ids)[j]
    y_j     = y[subject_ids == subject_id,] # observed y for subject j (vector)
    X_j     = x_mm[subject_ids == subject_id,] # fixed-effect design matrix for subject j (matrix)
    Z_j     = object$rand.covar[subject_ids == subject_id,]

    V_j <- Z_j %*% G %*% t(Z_j) + sigmae2 * diag(nrow(X_j))
    randeff <- (G %*% t(Z_j) %*% solve(V_j, y_j - X_j %*% beta))[,1]
    y_hat_random[subject_ids == subject_id,] <- Z_j %*% randeff
  }

  y_hat <- y_hat_fixed + y_hat_random
  pred_y <- NULL
  pred_y$y_hat <- y_hat
  pred_y$y_hat_random <- y_hat_random
  pred_y$y_hat_fixed <- y_hat_fixed
  class(pred_y) <- "predict.acml"
  pred_y
}


#' @exportS3Method
#' @importFrom stats residuals predict
residuals.acml <- function(object, digits=NULL, ...)
{
  y = object$response
  y_pred <- predict(object)
  resid_type1 <- y - y_pred$y_hat_fixed
  resid_type2 <- y - y_pred$y_hat
  resid <- NULL
  resid$resid_type1 <- resid_type1
  resid$resid_type2 <- resid_type2
  class(resid) <- "residuals.acml"
  resid
}

#' @exportS3Method
#' @importFrom graphics par
#' @importFrom grDevices rgb dev.flush dev.hold dev.interactive devAskNewPage
#' @importFrom stats qqline qqnorm rbinom
plot.acml <- function(
  x, digits=NULL, which=1:4,
  caption = list("Marginal Residuals",
                 "Marginal Residuals QQ",
                 "Conditional Residuals",
                 "Conditional Residuals QQ"),
  ask = prod(par("mfcol")) < length(which) && dev.interactive(),
  ...)
{
  object <- x
  y_pred <- predict(object)
  resid  <- residuals(object)
  oldpar <- par(ask = TRUE)
  on.exit(par(oldpar))

  show <- rep(FALSE, 4)
  show[which] <- TRUE

  if (ask)
  {
  	oask <- devAskNewPage(TRUE)
  	on.exit(devAskNewPage(oask))
  }

  # Type I residual plot
  if(show[1])
  {
    dev.hold()
    plot(y_pred$y_hat_fixed, resid$resid_type1,
         xlab = "Fitted (marginal)",
         ylab = "Residuals (Type I)",
         main = caption[1],
         ...)
    abline(h = 0, col = "blue")
    dev.flush()
  }

  if(show[2])
  {
    dev.hold()
    # Type I Q-Q plot
    qqnorm(resid$resid_type1, main = caption[2])
    qqline(resid$resid_type1, col = "blue")
    dev.flush()
  }

  # Type II residual plot
  if(show[3])
  {
    dev.hold()
    plot(y_pred$y_hat, resid$resid_type2,
         xlab = "Fitted (conditional)",
         ylab = "Residuals (Type II)",
         main = caption[3],
         ...)
    abline(h = 0, col = "red")
    dev.flush()
  }

  if(show[4])
  {
    dev.hold()
    # Type II Q-Q plot
    qqnorm(resid$resid_type2, main = caption[4])
    qqline(resid$resid_type2, col = "red")
    dev.flush()
  }
  invisible()
}

#' Retrieve the robust variance covariance matrix
#'
#' An S3 Method to retrieve the robust variance-covariance matrix of an ODS
#' model. It is derived via the sandwich estimator.
#'
#' @title robcov: Return the robust covariance matrix.
#' @param object to get robust variance covariance matrix.
#' @param ... other arguments
#' @export
robcov <- function(object, ...) UseMethod("robcov")

#' @exportS3Method
robcov.acml <- function(object, ...) object$robcov

#' @exportS3Method
logLik.acml <- function(object, ...)
  structure(object$LogL, nall = nrow(object$data),
    nobs = nrow(object$data), df=NA, class=c("logLik", "numeric"))

#' Fit model using ascertainment corrected likelihood model (ACML)
#'
#' Outcome dependent sampling designs need to be corrected when fitting a
#' statistical model for proper inferences. This routine will fit and return
#' the model fit using ascertainment corrected likelihood.
#'
#' @param formula `formula`; (or an object that can be
#'   coerced to that class): a symbolic description of the model to be fitted.
#'   The details of model specification are given under ‘Details’.
#' @param design an object of class 'odsdesign'. This specifies the
#'   design of sampling used for the fitting algorithm.
#' @param data `data.frame`; an optional data frame, list or environment (or
#'   object coercible by as.data.frame to a data frame) containing the variables
#'   in the model. If not found in data, the variables are taken from
#'   environment(formula), typically the environment from which acml is called.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process. (See additional details about how this
#'   argument interacts with data-dependent bases in the ‘Details’ section of
#'   the model.frame documentation.)
#' @param weights	an optional vector of weights to be used in the fitting
#'   process. Should be NULL or a numeric vector. If non-NULL, weighted least
#'   squares is used with weights weights (that is, minimizing sum(w*e^2));
#'   otherwise ordinary least squares is used. See also ‘Details’,
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param subset `logical`; an optional vector specifying a subset of
#'   observations to be used in the fitting process. (See additional details
#'   about how this argument interacts with data-dependent bases in the
#'   ‘Details’ section of the model.frame documentation.)
#' @param MI `logical(1)`; Is multiple imputation to be performed. Defaults to FALSE.
#' @param MImethod `character(1)`; Specifies multiple imputation method, defaults to 'direct'. Can also be 'indirect'.
#' @param na.action `function`; an optional vector specifying a subset of
#'   observations to be used in the fitting process. (See additional details
#'   about how this argument interacts with data-dependent bases in the
#'   ‘Details’ section of the model.frame documentation.)
#' @param verbose `numeric(1)`; Debugging information printing level
#'   (sent to the optimizer).
#' @param init `numeric`; Initial starting position for parameter search
#'   via the optimizer.
#' @param ... Optional additional parameters passed to sub methods.
#'
#' @export
#' @importFrom checkmate makeAssertCollection
#' @importFrom checkmate assert_formula assert_numeric assert_class assert_logical assert_choice
#' @importFrom checkmate reportAssertions
#' @importFrom stats lm model.matrix model.frame
#' @examples
#' data(gbti)
#' design <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
#'               data=gbti, quantiles=c(0.1, 0.9))
#' est <- acml(Response ~ Month*Genotype, design, gbti)
#' est
#' summary(est)
acml <- function(
  formula,
  design,
  data = NULL,
  subset = NULL,
  weights = NULL,
  MI = FALSE,
  MImethod = "direct",
  na.action = getOption('na.action'),
  verbose = 0L,
  init = NULL,
  ...)
{
  # Validate arguments
  coll <- makeAssertCollection()
  assert_formula(formula, add=coll)
  assert_class(design, "odsdesign", add=coll)
  assert_logical(MI, len=1, add=coll)
  assert_choice(MImethod, c("direct", "indirect"))
  reportAssertions(coll)

  # Duplicate of lm behavior
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  formula2 <- as.character(formula)
  formula2 <- as.formula(paste(formula[2], "~", formula[3], "+", design$id, collapse=''))
  mf[['formula']] <- formula2
  mf <- eval(mf, parent.frame())

  # Initial guess of coefficients
  if(is.null(init))
  {
    init <- c(coef(lm(formula, data=data, na.action=na.action)),
      c(1, 1, 1, 1)) # FIXME: Check for transformation
      # Intercept variance component, Slope Variance component, correlation, error variance
  }

  if(!is.numeric(mf[,design$id]) && !is.integer(mf[,design$id]))
    mf[,design$id] <- as.integer(as.factor(mf[,design$id]))

  if(!design$response %in% names(mf))
    stop("must have same response variable as design")
  if(!design$time %in% names(mf))
    stop("must have same time variable as design")

  mm <- model.matrix(formula2, mf, na.action=na.action)
  assert_true(all(as.character(mm[,design$id]) %in% names(design$p_sample_i)),
    .var.name='Group variables provided to acml that were not part of design', add=coll)
  reportAssertions(coll)

  fit <- acml_internal(
    y  = matrix(mf[,design$response]),
    x  = mm[,!(colnames(mm) %in% design$id)],
    z  = matrix(cbind(rep(1, nrow(mf)), mf[,design$time]), ncol=2),
    id = matrix(mf[,design$id]),
    w.function = if(design$method=='bivariate') 'bivar' else design$method,
    InitVals   = init,
    cutpoints  = as.vector(design$cutpoints),
    SampProb   = design$p_sample
  )
  fit$formula <- formula
  fit$design <- design
  fit$data   <- data
  fit$call   <- cl
  fit$model.matrix <- mm[,!(colnames(mm) %in% design$id)]
  fit$response     <- matrix(mf[,design$response])
  fit$rand.covar   <- matrix(cbind(rep(1, nrow(mf)), mf[,design$time]), ncol=2)
  fit$ids          <- matrix(mf[,design$id])

  class(fit) <- "acml"

  if(fit$Code == 3) warning("last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.")
  if(fit$Code == 4) warning("iteration limit exceeded.")
  if(fit$Code == 5) warning("maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.")

  fit
}
