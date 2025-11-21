library(MASS)
library(numDeriv)
library(knitr)
library(mvtnorm)
library(checkmate)
library(lme4)
library(tidyverse)

plot.odsdesign <- function(
    x,
    xlab   = "Intercept",
    ylab   = "Slope",
    main   = format(x$formula),
    sub    = paste(x$method, "design"),
    col    = if(x$method == 'mean') 'lightgray' else rgb(0,0,0,x$p_sample_i),
    lwd    = 1,
    lty    = 1,
    cutcol = 'red',
    cutlwd = 2,
    cutlty = 2,
    ...)
{
  if(x$method == 'mean')
  {
    if(xlab == 'Intercept') xlab <- "Mean"
    hist(x$z_i[1,], xlab=xlab, main=main, sub=sub,
         col=col, lwd=lwd, lty=lty, ...)
  } else
  {
    plot(x$z_i["intercept",], x$z_i["slope",],  # Lucy changed to rownames
         xlab=xlab, ylab=ylab,
         main=main, sub=sub,
         col=col, lwd=lwd, lty=lty, # Needed to prevent capture in call to lines below via ...
         ...)
  }

  if (x$method != 'bivariate')
  {
    for(i in colnames(x$cutpoints))
    {
      if(i %in% c('mean', 'intercept'))
        abline(v=x$cutpoints[,i], col=cutcol, lty=cutlty, lwd=cutlwd, ...) else
          abline(h=x$cutpoints[,i], col=cutcol, lty=cutlty, lwd=cutlwd, ...)
    }
  } else # bivariate
  {
    square <- function(x, y)
      lines(x[c(1, 1, 2, 2, 1)], y[c(1, 2, 2, 1, 1)],
            col=cutcol, lty=cutlty, lwd=cutlwd, ...)
    n <- nrow(x$cutpoints)
    for(i in 1:(n/2))
    {
      sel <- c(i, n+1-i)
      square(x$cutpoints[sel,1], x$cutpoints[sel,2])
    }
  }

  invisible(x)
}

transform_output <- function(InitVals, x = x, y = y, z = z){
  beta = InitVals[1:(ncol(x))] # note here are with beta_0, need to change to accommodate if no beta_0
  log_sigma_vc = InitVals[(ncol(x)+1):(ncol(x)+ncol(z))]
  log_inv_rho_vc = InitVals[(ncol(x)+ncol(z)+1):(ncol(x)+ncol(z) + choose(ncol(z), 2))]
  log_sigma_e = InitVals[(ncol(x)+ncol(z) + choose(ncol(z), 2)+1):length(InitVals)]
  sigma_vc = ifelse(is.infinite(exp(log_sigma_vc)), 100000, exp(log_sigma_vc))
  sigma_e = ifelse(is.infinite(exp(log_sigma_e)), 100000, exp(log_sigma_e))
  rho_vc = ifelse(is.infinite(exp(log_inv_rho_vc)), 1, (exp(log_inv_rho_vc) - 1)/(exp(log_inv_rho_vc) + 1))
  list(beta = beta, sigma_vc = sigma_vc, rho_vc = rho_vc, sigma_e = sigma_e)
}

#' @exportS3Method
model.matrix.odsdesign <- function(object, ...) as.matrix(object$model.frame)

#' @exportS3Method
print.odsdesign <- function(x, digits = max(3L, getOption("digits")), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), collapse="\n"),
      "\n\n",
      "Cutpoints:\n",
      sep="")
  print(round(x$cutpoints, digits=digits), ...)
  cat("\n")
  invisible(x)
}

#' @exportS3Method
#' @importFrom stats quantile
summary.odsdesign <- function(object, digits = max(3L, getOption("digits")), ...)
{
  ans <- object[c("call", "cutpoints")]
  ans$digits <- digits
  xx <- matrix(rep(NA, 24), ncol=4,dimnames=list(
    c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max." , "Mean"),
    c(names(object$model.frame)[1], rownames(object$z_i))
  ))
  xx[1:5,1] <- quantile(object$model.frame[,1], na.rm=TRUE, names=FALSE)
  xx[6,1]   <- mean(object$model.frame[,1], na.rm=TRUE)
  xx[1:5,2] <- quantile(object$z_i['mean',], na.rm=TRUE, names=FALSE)
  xx[6,2]   <- mean(object$z_i['mean',], na.rm=TRUE)
  xx[1:5,3] <- quantile(object$z_i['intercept',], na.rm=TRUE, names=FALSE)
  xx[6,3]   <- mean(object$z_i['intercept',], na.rm=TRUE)
  xx[1:5,4] <- quantile(object$z_i['slope',], na.rm=TRUE, names=FALSE)
  xx[6,4]   <- mean(object$z_i['slope',], na.rm=TRUE)

  ans$descriptive <- as.table(xx[c(1,2,6,3:5),])

  ans$N <- c(nrow(object$model.frame), ncol(object$z_i), sum(object$p_sample_i))
  names(ans$N) <- c("N", names(object$model.frame)[3], "E[N_sample]")

  class(ans) <- "summary.odsdesign"
  ans
}

#' @exportS3Method
print.summary.odsdesign <- function(x, digits = NULL, ...)
{
  if(is.null(digits)) digits <- x$digits
  cat("\nCall:\n",
      paste(deparse(x$call), collapse="\n"),
      "\n\n",
      "Cutpoints:\n",
      sep="")
  print(round(x$cutpoints, digits=digits), ...)
  cat("\n")
  print(round(x$descriptive, digits=digits), ...)
  cat("\n")
  print(round(x$N, digits=digits), ...)
  cat("\n")
  invisible(x)
}


ods <- function(
    formula,
    method,
    p_sample,
    data      = NULL,
    quantiles = NULL,
    cutpoints = NULL,
    subset    = NULL,
    weights   = NULL,
    na.action = getOption('na.action'),
    ProfileCol= NULL     ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
    # ... #FIXME: this doesn't work for now
)
{
  n_c <- length(p_sample)-1 # Number of cuts
  if(!is.null(method) && inherits(method, "character") && method == 'bivariate')
    n_c <- 2*n_c

  # Validate arguments
  coll <- makeAssertCollection()
  assert_formula(formula, add=coll)
  assert_numeric(p_sample,  add=coll, lower=0, upper=1, min.len=2, any.missing=FALSE)
  assert_numeric(quantiles, add=coll, lower=0, upper=1, len=length(p_sample)-1, null.ok=TRUE, any.missing=FALSE)
  assert_numeric(cutpoints, add=coll, len=n_c, null.ok=TRUE, any.missing=FALSE)
  assert_character(method,  add=coll, len=1, any.missing=FALSE, pattern="slope|intercept|bivariate|mean")
  assert_true(xor(is.null(quantiles), is.null(cutpoints)), .var.name="only one of quantiles or cutpoints can be specified", add=coll)
  assert_true(length(terms(formula))==3, .var.name="formula must have 3 terms", add=coll)
  assert_true(any(grepl("\\|", formula)), .var.name="formula must have an id specified, e.g. y~t|id", add=coll)
  reportAssertions(coll)

  # Square donut must have even number of quantiles or cutpoints
  if(method == 'bivariate')
  {
    assert_true(length(quantiles) %% 2 == 0, .var.name="length(quantiles) must be even for bivariate design", add=coll)
    assert_true(length(cutpoints) %% 2 == 0, .var.name="length(cutpoints) must be even for bivariate design", add=coll)
    reportAssertions(coll)
  }

  # Duplicate of lm behavior
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf[['formula']] <- as.formula(gsub("\\|", "+", format(formula)))
  mf <- eval(mf, parent.frame())

  if(!is.integer(mf[,3]) && !is.numeric(mf[,3]))
    mf[,3] <- as.numeric(as.factor(mf[,3]))

  # Construct z_i
  f <- if (ncol(mf) == 4) # Treat weights?
  {
    function(x) c(mean(x[,1]), coef(lm(x[,1]~x[,2], weights=x[,4], na.action=na.action)))
  } else
  {
    function(x) c(mean(x[,1]), coef(lm(x[,1]~x[,2], na.action=na.action)))
  }
  z_i <- sapply(split(mf, mf[,3]), f)
  rownames(z_i) <- c("mean", "intercept", "slope")

  z_mf <- model.matrix(reformulate(names(mf)[2], intercept = TRUE), data = mf)

  # Devise cutpoints if not specified
  if(is.null(cutpoints))
  {
    sel <- if(method == 'bivariate') c("intercept", "slope") else method
    cutpoints <- apply(z_i, 1, quantile, quantiles)[,sel, drop=FALSE]
  } else if(method == 'bivariate')
  {
    cutpoints <- matrix(
      cutpoints,
      nrow=2, byrow=FALSE,
      dimnames=list(1:2, c("intercept", "slope"))
    )
  } else
  {
    cutpoints <- matrix(
      cutpoints,
      nrow=2,
      dimnames=list(1:2, method))
  }

  # Create patient to probability
  p_sample_i <- if(method =='bivariate')
  {
    # Square donut(s)
    pmax(p_sample[as.numeric(
      cut(z_i['intercept',], c(-Inf, t(cutpoints[,'intercept']), Inf)))],
      p_sample[as.numeric(
        cut(z_i['slope',],     c(-Inf, t(cutpoints[,'slope']),     Inf)))])
  } else
  {
    p_sample[as.numeric(cut(z_i[method,], c(-Inf, t(cutpoints), Inf)))]
  }
  names(p_sample_i) <- colnames(z_i)

  smpl <- names(p_sample_i)[rbinom(length(p_sample_i), 1, p_sample_i) > 0]

  if(is.null(ProfileCol)){
    ProfileCol = NA
  }

  # Return design object
  structure(list(
    call        = cl,
    formula     = formula,
    model.frame = mf,
    method      = method,
    p_sample    = p_sample,
    p_sample_i  = p_sample_i,
    sample_ids  = smpl,
    response    = names(mf)[1],
    time        = names(mf)[2],
    id          = names(mf)[3],
    quantiles   = quantiles,
    cutpoints   = cutpoints,
    z_i         = z_i,
    z_mf        = z_mf,      #FIXME: only appropriate for one intercept and one slope
    weights     = weights,
    n_rand      = 2,     # Number of random effects, slope + intercept
    ProfileCol  = ProfileCol ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
  ),
  class="odsdesign"
  )
}


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

#' Calculate V_i = Z_i D t(Z_i) + sig_e^2 I_{n_i}
#'
#' Calculate V_i = Z_i D t(Z_i) + sig_e^2 I_{n_i}
#' @param zi n_i by q design matrix for the random effects
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return V_i
#' @export
#'
vi.calc <- function(zi, sigma.vc, rho.vc, sigma.e){
  SDMat.RE  <- diag(sigma.vc)
  ncolzi    <- ncol(zi) ## make sure this equals length(sigma.vc)
  nrowzi    <- nrow(zi)
  nERRsd    <- length(sigma.e)
  b         <- matrix(0,ncolzi,ncolzi)
  b[lower.tri(b, diag=FALSE)] <- rho.vc
  CorMat.RE <- t(b)+b+diag(rep(1,ncolzi))
  CovMat.RE <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
  zi %*% CovMat.RE %*% t(zi) + diag(rep(sigma.e^2, each=nrowzi/nERRsd))
}

#' Ascertainment correction piece for univariate sampling
#'
#' Calculate the (not yet log transformed) ascertainment correction under a univariate Q_i
#'
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 2)
#' @param SampProb Sampling probabilities from within each region (vector of length 3).
#' @param mu_q a scalar for the mean value of the Q_i distribution
#' @param sigma_q a scalar for the standard deviation of the Q_i distribution
#' @return Not yet log transformed ascertainment correction
#' @export
ACi1q <- function(cutpoints, SampProb, mu_q, sigma_q){
  CDFs <- pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
  sum( SampProb*(CDFs[2:length(CDFs)] - CDFs[1:(length(CDFs)-1)]) )
}

#' Ascertainment correction piece for bivariate sampling
#'
#' Calculate the (not yet log transformed) ascertainment correction under a bivariate Q_i
#'
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4: c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each of two sampling regions; central region and outlying region (vector of length 2).
#' @param mu_q a 2-vector for the mean value of the bivariate Q_i distribution.
#' @param sigma_q a 2 by 2 covariance matrix for the bivariate Q_i distribution.
#' @return Not yet log transformed ascertainment correction
#' @export
ACi2q <- function(cutpoints, SampProb, mu_q, sigma_q){
  (SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2]
}

#' Log of the Ascertainment correction for univariate sampling
#'
#' Calculate the log transformed ascertainment correction under a univariate Q_i.  Also return vi
#'
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param zi n_i by q design matric for random effects (intercept and slope)
#' @param wi the pre-multiplier of yi to generate the sampling variable q_i
#' @param beta mean model parameter vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 2)
#' @param SampProb Sampling probabilities from within each region (vector of length 3).
#' @return log transformed ascertainment correction
#' @export
logACi1q <- function(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb){
  vi      <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
  mu      <- xi %*% beta
  mu_q    <- (wi %*% mu)[,1]
  sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
  return(list(vi=vi, logACi=log(ACi1q(cutpoints, SampProb, mu_q, sigma_q))))
}


#' Log of the Ascertainment correction piece for bivariate sampling
#'
#' Calculate the log transformed ascertainment correction under a bivariate Q_i.  Also return vi
#'
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param zi n_i by q design matric for random effects (intercept and slope)
#' @param wi the pre-multiplier of yi to generate the sampling variable q_i
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).
#' @return log transformed ascertainment correction
#' @export
logACi2q <- function(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb){
  vi      <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
  mu      <- xi %*% beta
  mu_q    <- as.vector(wi %*% mu)
  sigma_q <- wi %*% vi %*% t(wi)
  ## for some reason the upper and lower triangles to not always equal, so I am taking their
  ## average.  Not sure this is a problem here
  ## but doing this to be safe.  Maybe can remove later once understood.
  sigma_q <- (sigma_q + t(sigma_q))/2
  #sigma_q[upper.tri(sigma_q)]  <- t(sigma_q)[upper.tri(sigma_q)]
  #sigma_q[2,1] <- sigma_q[1,2]
  return(list(vi=vi, logACi= log( ACi2q(cutpoints=cutpoints, SampProb=SampProb, mu_q=mu_q, sigma_q=sigma_q))))
}


#' Calculate a subject-specific contribution to a log-likelihood for longitudinal normal data
#'
#' Calculate a subject-specific contribution to a log-likelihood for longitudinal normal data
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param beta mean model parameter vector
#' @param vi the variance covariance matrix (ZDZ+Sige2*I)
#' @return subject specific contribution to the log-likelihood
#' @export
#'
li.lme <- function(yi, xi, beta, vi){
  resid <- yi - xi %*% beta
  -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}

#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by q design matrix for random effects
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivariate."  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param Keep.liC If FALSE, the function returns the conditional log likelihood across all subjects.  If TRUE, subject specific contributions and exponentiated subject specific ascertainment corrections are returned in a list.
#' @return If Keep.liC=FALSE, conditional log likelihood.  If Keep.liC=TRUE, a two-element list that contains subject specific likelihood contributions and exponentiated ascertainment corrections.
#' @export
#'
LogLikeC2 <- function(y, x, z, w.function, id, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, Weights, Keep.liC=FALSE){

  subjectData    <- CreateSubjectData(id=id,y=y,x=x,z=z,Weights=Weights,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function)
  liC.and.logACi <- lapply(subjectData, LogLikeiC2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)

  if (Keep.liC == FALSE){out <- -1*Reduce('+', liC.and.logACi)[1]  ## sum ss contributions to liC
  }else{ out <- list(liC    = c(unlist(sapply(liC.and.logACi, function(x) x[1]))), ## ss contributions to liC
                     logACi = c(unlist(sapply(liC.and.logACi, function(x) x[2]))))} ## ss contributions to ACi
  out
}



#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i
#' @param w.function options include "mean" "intercept" "slope" and "bivariate"
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).
#' @return ss contributions to the conditional log likelihood.  This is an internal function used by LogLikeC2
#' @export
#'
#'
LogLikeiC2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){
  yi          <- subjectData[["yi"]]
  xi          <- subjectData[["xi"]]
  zi          <- subjectData[["zi"]]
  Weights.i   <- subjectData[["Weights.i"]]
  w.function  <- subjectData[["w.function.i"]]
  SampProb    <- subjectData[["SampProb.i"]]
  cutpoints   <- subjectData[["cutpoints.i"]]
  ni          <- length(yi)
  t.zi        <- t(zi)
  ##########
  ##########
  ##########
  ##########
  wi.tmp <- solve(t.zi %*% zi) %*% t.zi
  if (!(w.function %in% c("bivariate", "mvints", "mvslps"))){
    if (w.function %in% c("intercept", "intercept1")){ wi <- wi.tmp[1,]
    } else if (w.function %in% c("slope", "slope1")){  wi <- wi.tmp[2,]
    } else if (w.function %in% c("intercept2")){       wi <- wi.tmp[3,]
    } else if (w.function %in% c("slope2")){           wi <- wi.tmp[4,]
    } else if (w.function=="mean"){                    wi <- t(rep(1/ni, ni))
    }
    wi         <- matrix(wi, 1, ni)
    tmp        <- logACi1q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
    logACi     <- tmp[["logACi"]]
    liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  }else {
    if (w.function %in% c("bivariate")){ wi <- wi.tmp[c(1,2),]
    } else if (w.function %in% c("mvints")){ wi <- wi.tmp[c(1,3),]
    } else if (w.function %in% c("mvslps")){ wi <- wi.tmp[c(2,4),]
    }
    tmp        <- logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
    logACi     <- tmp[["logACi"]]
    liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  }
  ##########
  ##########
  ##########
  ##########
  # if (w.function != "bivariate"){
  #     if (w.function %in% c("intercept", "intercept1")){ wi<- (solve(t.zi %*% zi) %*% t.zi)[1,]
  #     } else if (w.function %in% c("slope", "slope1")){     wi<- (solve(t.zi %*% zi) %*% t.zi)[2,]
  #     } else if (w.function %in% c("intercept2")){ wi<- (solve(t.zi %*% zi) %*% t.zi)[3,]
  #     } else if (w.function %in% c("slope2")){     wi<- (solve(t.zi %*% zi) %*% t.zi)[4,]
  #     } else if (w.function=="mean"){     wi <- t(rep(1/ni, ni))
  #     }
  #     wi         <- matrix(wi, 1, ni)
  #     tmp        <- logACi1q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
  #     logACi     <- tmp[["logACi"]]
  #     liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  #
  # }else{
  #     wi         <- solve(t.zi %*% zi) %*% t.zi
  #     tmp        <- logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
  #     logACi     <- tmp[["logACi"]]
  #     liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  #
  # }
  return(c(liC, logACi))
}


#' Gradient of the log of the ascertainment correction piece for sampling based on bivariate Q_i.
#'
#' Calculate the gradient of the log transformed ascertainment correction under designs that sample based on a bivariate Q_i (numerically).
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function, SampProb, cutpoints
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param sigmae std dev of the measurement error distribution
#' @return gradient of the log transformed ascertainment correction under the bivariate sampling design
#' @importFrom numDeriv grad
#' @importFrom mvtnorm pmvnorm
#' @export
logACi2q.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){

  yi          <- subjectData[["yi"]]
  xi          <- subjectData[["xi"]]
  zi          <- subjectData[["zi"]]
  Weights.i   <- subjectData[["Weights.i"]]  ## not used
  w.function  <- subjectData[["w.function.i"]]
  SampProb    <- subjectData[["SampProb.i"]]
  cutpoints   <- subjectData[["cutpoints.i"]]

  t.zi        <- t(zi)
  #wi          <- solve(t.zi %*% zi) %*% t.zi
  ###########################
  ###########################
  ###########################
  ###########################
  wi.tmp      <- solve(t.zi %*% zi) %*% t.zi
  if (w.function %in% c("bivariate")){         wi <- wi.tmp[c(1,2),]
  } else if (w.function %in% c("mvints")){ wi <- wi.tmp[c(1,3),]
  } else if (w.function %in% c("mvslps")){ wi <- wi.tmp[c(2,4),]
  }
  ###########################
  ###########################
  ###########################
  ###########################
  t.wi        <- t(wi)

  param   <- c(beta, sigma.vc, rho.vc, sigma.e)
  npar    <- length(param)

  len.beta     <- length(beta)
  len.sigma.vc <- length(sigma.vc)
  len.rho.vc   <- length(rho.vc)
  len.sigma.e  <- length(sigma.e)

  beta.index   <- c(1:len.beta)
  vc.sd.index  <- len.beta + (c(1:len.sigma.vc))
  vc.rho.index <-  len.beta + len.sigma.vc + (c(1:len.rho.vc))
  err.sd.index <- len.beta + len.sigma.vc + len.rho.vc + c(1:len.sigma.e)

  Deriv   <- sapply(1:npar,  function(rr)
  {
    grad(function(x) { new.param <- param
    new.param[rr] <- x
    vi      <- vi.calc(zi, new.param[vc.sd.index], new.param[vc.rho.index], new.param[err.sd.index])
    mu_q    <- as.vector(wi %*% (xi %*% new.param[c(beta.index)]))
    sigma_q <- wi %*% vi %*% t.wi
    ## for some reason the upper and lower triangles to not always equal.  Not sure this is a problem here
    ## but doing this to be safe.  Maybe can remove later once understood.
    sigma_q <- (sigma_q + t(sigma_q))/2
    pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]]
    },
    param[rr],
    method="simple",
    method.args=list(eps=1e-4))
  }
  )

  vi      <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
  mu_q    <- as.vector(wi %*% (xi %*% beta))
  sigma_q <- wi %*% vi %*% t.wi
  ## for some reason the upper and lower triangles to not always equal, so I am forcing
  ## the upper triangle to equal the lower triangle.
  sigma_q <- (sigma_q + t(sigma_q))/2

  (SampProb[1]-SampProb[2])*Deriv / ACi2q(cutpoints, SampProb, mu_q, sigma_q)
}



#' Gradient of the log of the ascertainment correction piece for sampling based on univariate Q_i
#'
#' Calculate the gradient of the log transformed ascertainment correction for sampling based on univariate Q_i
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function.i, SampProb.i, cutpoints.i
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return gradient of the log transformed ascertainment correction under univariate $Q_i$
#' @export
#' @importFrom stats dnorm
logACi1q.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){

  yi          <- subjectData[["yi"]]
  xi          <- subjectData[["xi"]]
  zi          <- subjectData[["zi"]]
  Weights.i   <- subjectData[["Weights.i"]]  ## not used
  w.function  <- subjectData[["w.function.i"]]
  SampProb    <- subjectData[["SampProb.i"]]
  cutpoints   <- subjectData[["cutpoints.i"]]

  t.zi        <- t(zi)
  ni          <- length(yi)
  #if (w.function=="mean")      wi <- t(rep(1/ni, ni))
  #if (w.function=="intercept") wi<- (solve(t.zi %*% zi) %*% t.zi)[1,]
  #if (w.function=="slope")     wi<- (solve(t.zi %*% zi) %*% t.zi)[2,]
  if (w.function %in% c("intercept", "intercept1")){ wi<- (solve(t.zi %*% zi) %*% t.zi)[1,]
  } else if (w.function %in% c("slope", "slope1")){     wi<- (solve(t.zi %*% zi) %*% t.zi)[2,]
  } else if (w.function %in% c("intercept2")){ wi<- (solve(t.zi %*% zi) %*% t.zi)[3,]
  } else if (w.function %in% c("slope2")){     wi<- (solve(t.zi %*% zi) %*% t.zi)[4,]
  } else if (w.function=="mean"){     wi <- t(rep(1/ni, ni))
  }

  wi      <- matrix(wi, 1, ni)

  vi        <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
  t.wi      <- t(wi)
  wi.zi     <- wi %*% zi
  t.wi.zi   <- t(wi.zi)

  mu        <- xi %*% beta
  mu_q      <- (wi %*% mu)[,1]
  sigma_q   <- sqrt((wi %*% vi %*% t.wi)[1,1])

  l <- ACi1q(cutpoints, SampProb, mu_q, sigma_q)
  p <- SampProb[1:(length(SampProb)-1)] - SampProb[2:(length(SampProb))]
  f <- dnorm(cutpoints, mu_q, sigma_q)

  d_li_beta <- (wi %*% xi) * sum(p*f) / l
  f_alpha_k <- sum(p*f*(cutpoints - mu_q)) / (l * 2* sigma_q^2 )

  #################################################################################
  ### now calculate d_li_dalpha
  #################################################################################
  len.sigma.vc <- length(sigma.vc)
  len.rho.vc   <- length(rho.vc)
  len.sigma.e  <- length(sigma.e)

  SDMat.RE     <- diag(sigma.vc)

  b0         <- matrix(0, len.sigma.vc, len.sigma.vc)
  b0[lower.tri(b0, diag=FALSE)] <- rho.vc
  tbb0       <- t(b0) + b0
  CorMat.RE  <- tbb0 +  diag(rep(1,len.sigma.vc))

  D             <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
  m             <- matrix(0,len.sigma.vc, len.sigma.vc)

  ## derivatives w.r.t variance components SDs
  dVi.dsigma.vc <- NULL
  for (mmm in 1:len.sigma.vc){
    m1               <- m
    m1[mmm,]         <- m1[,mmm] <- 1
    m1[mmm,mmm]      <- 2
    tmp              <- sigma.vc[mmm]+ 1*(sigma.vc[mmm]==0) ## to prevent unlikely division by 0
    dViMat.dsigma.vc <- m1 * D / tmp
    dVi.dsigma.vc    <- c(dVi.dsigma.vc, (wi.zi %*% dViMat.dsigma.vc %*% t.wi.zi)[1,1])
  }

  ## derivatives w.r.t variance components rhos
  b1         <- matrix(0,len.sigma.vc,len.sigma.vc)
  b1[lower.tri(b1, diag=FALSE)] <-  c(1:len.rho.vc)
  tbb1       <- t(b1) + b1

  dVi.drho.vc <- NULL
  for (mmm in 1:len.rho.vc){
    tmp                   <- which(tbb1==mmm, arr.ind=TRUE)
    m2                    <- m
    m2[tmp[1,1],tmp[1,2]] <- 1
    m2[tmp[2,1],tmp[2,2]] <- 1
    tmp                   <- rho.vc[mmm] + 1*(rho.vc[mmm]==0) ## to prevent division by 0
    dViMat.drho.vc        <- m2 * D / tmp
    dVi.drho.vc           <- c(dVi.drho.vc, (wi.zi %*% dViMat.drho.vc %*% t.wi.zi)[1,1])
  }

  ## derivatives w.r.t error sds
  dVi.dsigma.e <- NULL
  for (mmm in 1:len.sigma.e){
    dsigma.e.vec       <- 2*sigma.e
    dsigma.e.vec[-mmm] <- 0
    dViMat.dsigma.e    <- diag(rep(dsigma.e.vec, each=ni/len.sigma.e))
    dVi.dsigma.e       <- c(dVi.dsigma.e, (wi %*% dViMat.dsigma.e %*% t.wi)[1,1])
  }

  c(d_li_beta, c(f_alpha_k * c(dVi.dsigma.vc, dVi.drho.vc, dVi.dsigma.e)))
}


#' Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))
#'
#' Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))
#'
#' @param subjectData a list that contains yi, xi, zi, Weights.i.  Note that Weights.i is used for inverse probability weighting only.
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return Subject specific contribution to the log-likelihood score (also returns marginal Vi=Cov(Y|X))
#' @export
li.lme.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){
  yi          <- subjectData[["yi"]]
  xi          <- subjectData[["xi"]]
  zi          <- subjectData[["zi"]]
  t.zi        <- t(zi)
  ni          <- length(yi)
  Weights.i   <- subjectData[["Weights.i"]]

  resid         <- yi - xi %*% beta
  vi            <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
  inv.v         <- solve(vi)
  t.resid       <- t(resid)
  t.resid.inv.v <- t.resid %*% inv.v
  inv.v.resid   <- inv.v %*% resid

  dli.dbeta <- t(xi) %*% inv.v.resid # for Beta

  len.sigma.vc <- length(sigma.vc)
  len.rho.vc   <- length(rho.vc)
  len.sigma.e  <- length(sigma.e)

  SDMat.RE     <- diag(sigma.vc)

  b0         <- matrix(0, len.sigma.vc, len.sigma.vc)
  b0[lower.tri(b0, diag=FALSE)] <- rho.vc
  tbb0       <- t(b0) + b0
  CorMat.RE  <- tbb0 +  diag(rep(1,len.sigma.vc))

  D             <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
  m             <- matrix(0,len.sigma.vc, len.sigma.vc)

  ## derivatives w.r.t variance components SDs
  dli.dsigma.vc <- NULL
  for (mmm in 1:len.sigma.vc){
    m1               <- m
    m1[mmm,]         <- m1[,mmm] <- 1
    m1[mmm,mmm]      <- 2
    tmp              <- sigma.vc[mmm]+ 1*(sigma.vc[mmm]==0) ## to prevent unlikely division by 0
    dViMat.dsigma.vc <- m1 * D / tmp
    tmp              <- zi %*% dViMat.dsigma.vc %*% t.zi
    dli.dsigma.vc    <- c(dli.dsigma.vc, -0.5*(sum(diag(inv.v %*% tmp)) - t.resid.inv.v %*% tmp %*% inv.v.resid)[1,1])
  }

  ## derivatives w.r.t variance components rhos
  b1         <- matrix(0,len.sigma.vc,len.sigma.vc)
  b1[lower.tri(b1, diag=FALSE)] <-  c(1:len.rho.vc)
  tbb1       <- t(b1) + b1

  dli.drho.vc <- NULL
  for (mmm in 1:len.rho.vc){
    tmp                   <- which(tbb1==mmm, arr.ind=TRUE)
    m2                    <- m
    m2[tmp[1,1],tmp[1,2]] <- 1
    m2[tmp[2,1],tmp[2,2]] <- 1
    tmp                   <- rho.vc[mmm] + 1*(rho.vc[mmm]==0) ## to prevent division by 0
    dViMat.drho.vc        <- m2 * D / tmp
    tmp                   <- zi %*% dViMat.drho.vc %*% t.zi
    dli.drho.vc           <- c(dli.drho.vc, -0.5*(sum(diag(inv.v %*% tmp)) - t.resid.inv.v %*% tmp %*% inv.v.resid)[1,1])
  }

  ## derivatives w.r.t error sds
  dli.dsigma.e <- NULL
  for (mmm in 1:len.sigma.e){
    dsigma.e.vec       <- 2*sigma.e
    dsigma.e.vec[-mmm] <- 0
    tmp                <- diag(rep(dsigma.e.vec, each=ni/len.sigma.e))
    dli.dsigma.e       <- c(dli.dsigma.e, -0.5*(sum(diag(inv.v %*% tmp)) - t.resid.inv.v %*% tmp %*% inv.v.resid)[1,1])
  }
  #print(c(dli.dsigma.vc, dli.drho.vc, dli.dsigma.e))
  list(gr=append(dli.dbeta, c(dli.dsigma.vc, dli.drho.vc, dli.dsigma.e))*Weights.i,
       vi=vi)
}


#' Calculate the gradient of the conditional likelihood
#' @description Calculate the gradient of the conditional likelihood for the univariate and bivariate sampling cases across all subjects (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if CheeseCalc=TRUE.
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[1,])}, "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[1,])}. "intercept2" (second intercept of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[3,])}, "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[4,])} "bivariate" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[c(1,2),])} "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])} "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])}.  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param CheeseCalc If FALSE, the function returns the gradient of the conditional log likelihood across all subjects.  If TRUE, the cheese part of the sandwich esitmator is calculated.
#' @return If CheeseCalc=FALSE, gradient of conditional log likelihood.  If CheeseCalc=TRUE, the cheese part of the sandwich estimator is calculated.
#' @export
LogLikeC.Score2 <- function(y, x, z, w.function, id, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, Weights, CheeseCalc=FALSE){
  param.vec <- c(beta, log(sigma.vc),log((1+rho.vc)/(1-rho.vc)),log(sigma.e))
  #print(c("blahblah", param.vec))
  npar     <- length(param.vec)

  len.beta     <- length(beta)
  len.sigma.vc <- length(sigma.vc)
  len.rho.vc   <- length(rho.vc)
  len.sigma.e  <- length(sigma.e)

  beta.index    <- c(1:len.beta)
  vc.sd.index   <- len.beta + (c(1:len.sigma.vc))
  vc.rho.index  <-  len.beta + len.sigma.vc + (c(1:len.rho.vc))
  err.sd.index  <- len.beta + len.sigma.vc + len.rho.vc + c(1:len.sigma.e)
  notbeta.index <- c(vc.sd.index,vc.rho.index,err.sd.index)

  subjectData <- CreateSubjectData(id=id,y=y,x=x,z=z,Weights=Weights,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function)

  UncorrectedScorei <- lapply(subjectData, li.lme.score2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)
  Gradienti         <- lapply(UncorrectedScorei, function(x) x[['gr']]) ## create a list of ss contributions to gradient
  UncorrectedScore  <- Reduce('+', Gradienti)  ## Note if using IPW this is actually a corrected score (corrected by the IPW)
  #print(c("blah1", UncorrectedScore))
  ## NOTE HERE: I used the first element of w.function in this call.  This means, for now, we cannot mix bivariate with other
  ## sampling schemes.  This also applies to the cheese calculation
  if (!(w.function[[1]] %in% c("bivariate","mvints","mvslps"))){
    logACi.Score <- lapply(subjectData, logACi1q.score2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)
    logAC.Score  <- Reduce('+', logACi.Score)
    CorrectedScore <- UncorrectedScore + logAC.Score
  }else{
    logACi.Score <- lapply(subjectData, logACi2q.score2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)
    logAC.Score  <- Reduce('+', logACi.Score)
    CorrectedScore <- UncorrectedScore - logAC.Score  ## notice this has the opposite sign compared to above.  Remember to check
  }
  #print(c("blah2", CorrectedScore))
  if (CheeseCalc==TRUE){
    if (!(w.function[[1]] %in% c("bivariate","mvints","mvslps"))){ GradiMat <- mapply("+", Gradienti, logACi.Score)
    }else{                           GradiMat <- mapply("-", Gradienti, logACi.Score)} ## notice this has the opposite sign compared to above.  Remember to check
    ## Need to use the chain rule: note that param,vec is on the unconstrained scale but Gradi was calculated on the constrained parameters
    GradiMat[notbeta.index,] <- GradiMat[notbeta.index,]*c(exp(param.vec[vc.sd.index]), 2*exp(param.vec[vc.rho.index])/((exp(param.vec[vc.rho.index])+1)^2), exp(param.vec[err.sd.index]))
    cheese <- matrix(0,  npar, npar)
    for (mm in 1:ncol(GradiMat)) cheese <- cheese + outer(GradiMat[,mm], GradiMat[,mm])
  }
  out <- -CorrectedScore
  if (CheeseCalc==TRUE) out <- cheese
  out
}

#' Calculate the ascertainment corrected log likelihood and score
#' @description
#' Calculate the ascertainment corrected log likelihood and score for different designs
#' @param params parameter vector c(beta, log(sigma0), log(sigma1), rho, sigmae)
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param id sum(n_i) vector of subject ids
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[1,])}, "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[1,])}. "intercept2" (second intercept of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[3,])}, "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[4,])} "bivariate" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \eqn{(solve(t.zi %*% zi) %*% t.zi)[c(1,2),])} "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])} "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \eqn{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])}.  There should be one unique value per subject
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maimizing the log likelihood for all other parameters while fixing these columns at the values of params at the location of ProfileCol
#' @param Keep.liC If TRUE outputs subject specific conditional log lileihoods to be used for the imputation procedure described in the AOAS paper keep z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @return The conditional log likelihood with a "gradient" attribute (if Keep.liC=FALSE) and subject specific contributions to the conditional likelihood if Keep.liC=TRUE).
#' @export
LogLikeCAndScore2 <- function(params, y, x, z, id, w.function, cutpoints, SampProb, Weights, ProfileCol=NA, Keep.liC=FALSE){
  npar   <- length(params)

  nbeta <- ncol(x)
  nVCsd <- ncol(z)
  nVCrho <- choose(nVCsd,2)
  nERRsd <- npar-nbeta-nVCsd-nVCrho

  beta.index   <- c(1:nbeta)
  vc.sd.index  <- nbeta + (c(1:nVCsd))
  vc.rho.index <- nbeta + nVCsd + (c(1:nVCrho))
  err.sd.index <- nbeta + nVCsd + nVCrho + c(1:nERRsd)

  beta     <- params[beta.index]
  sigma.vc <- exp(params[vc.sd.index])
  rho.vc   <- (exp(params[vc.rho.index])-1) / (exp(params[vc.rho.index])+1)
  sigma.e  <- exp(params[err.sd.index])

  out     <- LogLikeC2( y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e, cutpoints=cutpoints,
                        SampProb=SampProb, Weights=Weights, Keep.liC=Keep.liC)
  GRAD    <- LogLikeC.Score2(y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e, cutpoints=cutpoints,
                             SampProb=SampProb, Weights=Weights)
  ## Need to use the chain rule: note that params is on the unconstrained
  ## scale but GRAD was calculated on the constrained parameters
  GRAD[vc.sd.index]  <- GRAD[vc.sd.index]*exp(params[vc.sd.index])
  GRAD[vc.rho.index] <- GRAD[vc.rho.index]*2*exp(params[vc.rho.index])/((exp(params[vc.rho.index])+1)^2)
  GRAD[err.sd.index] <- GRAD[err.sd.index]*exp(params[err.sd.index])
  ## Force the gradient of the fixed parameter to be zero, so that it does not move
  if (!is.na(ProfileCol)) GRAD[ProfileCol] <- 0
  attr(out,"gradient") <- GRAD

  out
}

#' Create a list of subject-specific data
#'
#' @param id sum(n_i) vector of subject ids
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivariate."  There should be one unique value per subject
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @export
CreateSubjectData <- function(id,y,x,z,Weights,SampProb,cutpoints,w.function){
  id.tmp        <- split(id,id)
  y.tmp         <- split(y,id)
  x.tmp         <- split(x,id)
  z.tmp         <- split(z,id)
  Weights.tmp <- split(Weights,id)
  SampProb.tmp  <- split(SampProb,id)
  cutpoints.tmp  <- split(cutpoints,id)
  w.function.tmp  <- split(w.function,id)

  subjectData <- vector('list', length=length(unique(id)))
  subjectData <- list()
  uid <- as.character(unique(id))
  for(j in seq(along=uid)){
    i <- uid[j]
    subjectData[[j]] <- list(idi          = as.character(unique(id.tmp[[i]])),
                             xi           = matrix(x.tmp[[i]], ncol=ncol(x)),
                             zi           = matrix(z.tmp[[i]], ncol=ncol(z)),
                             yi           = y.tmp[[i]],
                             Weights.i    = unique(Weights.tmp[[i]]),
                             SampProb.i   = matrix(SampProb.tmp[[i]], ncol=ncol(SampProb))[1,],
                             w.function.i = as.character(unique(w.function.tmp[[i]])),
                             cutpoints.i  = matrix(cutpoints.tmp[[i]], ncol=ncol(cutpoints))[1,])
  }
  names(subjectData) <- uid
  subjectData
}

## If you do not want to use the ascertainment correction term in the conditional likelihood
## set all SampProb values equal to each other.  This would be the case if you were doing
## straightforward maximum likelihood (albeit computationally inefficient) or weighted likelihood.

#' Fitting function: ACML or WL for a linear mixed effects model (random intercept and slope)
#'
#' @param formula.fixed formula for the fixed effects (of the form y~x)
#' @param formula.random formula for the random effects (of the form ~z).  Right now this model only fits random intercept and slope models.
#' @param data data frame that should contain everything in formula.fixed, formula.random, id, and Weights.  It does not include: w.function, cutpoints, SampProb
#' @param id sum(n_i) vector of subject ids (a variable contained in data)
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi* zi) * t.zi)[1,]), "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)[1,]). "intercept2" (second intercept of the regression of the Yi ~
##zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)[3,]), "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)[2,]), "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)[2,]), "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)[4,]) "bivariate" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)[c(1,2),]) "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) * t.zi)[c(1,3),]) "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)[c(1,3),]).  There should be one unique value per subject.  NOTE: We also have the same designs but for BLUP based sampling, in which case, the character string should begin with "blup.".  For example "blup.intercept". There should be one unique value per subject but n_i replicates of that value.  Note that w.function should NOT be in the dat dataframe.
#' @param InitVals starting values for c(beta, log(sigma0), log(sigma1), log((1+rho)/(1-rho)), log(sigmae))
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions for individual subjects.  If using a low, medium, high, sampling scheme, this is a sum(n_i) by 2 matrix that must be a distinct object not contained in the dat dataframe.  Each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high.  If using a square doughnut design this should be sum(n_i) by 4 matrix (var1lower, var1upper, var2lower, var2upper). Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region. For low medium high sampling, each row is a vector of length 3 with sampling probabilities for each region. For bivariate stratum sampling each row is a vector of length 2 with sampling probabilities for the inner and outer strata. Each subject should have n_i rows of the same values.  Not in data.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  This should be a variable in the data dataframe. It should only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar.  If not doing IPWL, this must be a vectors of 1s.
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maimizing the log likelihood for all other parameters
#'                   while fixing these columns at the values of InitVals[ProfileCol]
#' @return Ascertainment corrected Maximum likelihood: Ests, covar, logLik, code, robcov
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats na.omit
#' @importFrom stats nlm
#' @importFrom stats pnorm
#' @export

acml_internal <- function(formula,
                          design,
                          data,
                          InitVals
){ ## only used for blup sampling
  if(is.null(formula)) {stop('Specify the formula of the model.  It is currently NULL.')}
  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
    warning('data converted to data.frame.')
  }
  terms = unique( c(all.vars(formula),design$id,design$weights) )
  data  = data[,terms]

  if(any(is.na(data))) data = na.omit(data)

  id    = data$id = data[ , design$id]

  fixed.f = model.frame(formula, data)
  fixed.t = attr(fixed.f, "terms")
  y       = model.response(fixed.f,'numeric')
  uy      = unique(y)
  x       = model.matrix(formula, fixed.f)
  z_id    = data.frame(id = design$model.frame[,design$id], design$z_mf)
  z       = as.matrix(z_id[z_id$id %in% unique(id),-1])

  #if (is.na(SampProb[1])) SampProb = c(1,1,1)
  if (is.null(design$weights)){
    Weights    = data$Weights = rep(1, nrow(data))    #for non-IPW case
  }else{
    Weights0   = design$weights
    Weights    = data$Weights = data[ , Weights0]
  }

  if (design$method == "bivariate"){
    cutpoints   <- matrix(rep(design$cutpoints, length(y)), ncol = 4, byrow = T)
    SampProb    <- matrix(rep(design$p_sample[1:2], length(y)), ncol = 2, byrow = T)
  } else {
    cutpoints   <- matrix(rep(design$cutpoints, length(y)), ncol = 2, byrow = T)
    SampProb    <- matrix(rep(design$p_sample, length(y)), ncol = 3, byrow = T)
  }
  w.function  <- rep(design$method, length(y))

  acml.fit <- nlm(LogLikeCAndScore2,
                  InitVals,
                  y=y,
                  x=x,
                  z=z,
                  id=id,
                  w.function=w.function,
                  cutpoints=cutpoints,
                  SampProb=SampProb,
                  Weights=Weights,
                  ProfileCol=design$ProfileCol,
                  stepmax=4, iterlim=250,
                  check.analyticals = TRUE, print.level=0)

  ## Calculate the observed information and then invert to get the covariance matrix
  npar        <- length(acml.fit$estimate)
  Hessian.eps <- 1e-7
  eps.mtx     <- diag(rep(Hessian.eps, npar))
  grad.at.max <- acml.fit$gradient
  ObsInfo.tmp <- ObsInfo <- matrix(NA, npar, npar)


  ## Observed Information## Observed Informatiyon
  for (j in 1:npar){
    temp            <- LogLikeCAndScore2(acml.fit$estimate+eps.mtx[j,],
                                         y=y,
                                         x=x,
                                         z=z,
                                         id=id,
                                         w.function=w.function,
                                         cutpoints=cutpoints,
                                         SampProb=SampProb,
                                         Weights=Weights,
                                         ProfileCol=design$ProfileCol)
    ObsInfo.tmp[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
  }
  for (m in 1:npar){
    for (n in 1:npar){ ObsInfo[m,n] <-  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}

  ## Cheese part of the sandwich estimator
  nbeta <- ncol(x)
  nVCsd <- ncol(z)
  nVCrho <- choose(nVCsd,2)
  nERRsd <- npar-nbeta-nVCsd-nVCrho

  beta.index   <- c(1:nbeta)
  vc.sd.index  <- nbeta + (c(1:nVCsd))
  vc.rho.index <- nbeta + nVCsd + (c(1:nVCrho))
  err.sd.index <- nbeta + nVCsd + nVCrho + c(1:nERRsd)

  Cheese <- LogLikeC.Score2(y=y,
                            x=x,
                            z=z,
                            w.function=w.function,
                            cutpoints=cutpoints,
                            SampProb=SampProb,
                            id=id, beta=acml.fit$estimate[beta.index],
                            sigma.vc=exp(acml.fit$estimate[vc.sd.index]),
                            rho.vc=   (exp(acml.fit$estimate[vc.rho.index])-1) / (exp(acml.fit$estimate[vc.rho.index])+1),
                            sigma.e=exp(acml.fit$estimate[err.sd.index]),
                            Weights=Weights,
                            CheeseCalc=TRUE)

  if (!is.na(design$ProfileCol)){
    acml.fit$estimate <- acml.fit$estimate[-design$ProfileCol]
    ObsInfo           <- ObsInfo[-ProfileCol, -design$ProfileCol]
    Cheese            <- Cheese[-ProfileCol, -design$ProfileCol]
  }


  out              <- NULL
  out$call         <- match.call()
  out$coefficients <- acml.fit$estimate
  out$covariance   <- solve(ObsInfo)
  out$robcov       <- solve(ObsInfo)%*%Cheese%*%solve(ObsInfo)
  out$logLik       <- -acml.fit$minimum
  out$Code         <- acml.fit$code
  attr(out,'args') <- list(formula    = formula,
                           design_formula = design$formula,
                           id         = id,
                           w.function = w.function,
                           cutpoints  = cutpoints,
                           SampProb   = SampProb,
                           Weights    = Weights,
                           WeightsVar = design$weights,
                           ProfileCol = design$ProfileCol)
  if(kappa(out$covar) > 1e5) warning("Poorly Conditioned Model")
  out
}


# S3 methods for acml

triangle <- function(n) n*(n+1)/2

#' @export
#' @rdname coef
fixef <- function(object, ...) UseMethod("fixef")

#' @export
fixef.acml <- function(object, ...)
{
  est <- object$coefficients
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
  est <- object$coefficients
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
    y_j     = y[subject_ids == subject_id] # observed y for subject j (vector)
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
  structure(object$logLik, nall = nrow(object$data),
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
    data,
    subset = NULL,
    # weights = NULL,
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

  if(!is.numeric(data[,design$id]) && !is.integer(data[,design$id]))
    data[,design$id] <- as.integer(as.factor(data[,design$id]))

  if(!design$response %in% names(data))
    stop("must have same response variable as design")
  if(!design$time %in% names(data))
    stop("must have same time variable as design")

  mm <- model.matrix(formula, data, na.action=na.action)

  # mm <- model.matrix(formula2, mf, na.action=na.action)
  assert_true(all(as.character(data[,design$id]) %in% names(design$p_sample_i)),
              .var.name='Group variables provided to acml that were not part of design', add=coll)
  reportAssertions(coll)


  fit <- acml_internal(
    formula  = formula,
    design   = design,
    InitVals = init,
    data     = data
  )
  fit$formula <- formula
  fit$design <- design
  fit$data   <- data
  fit$call   <- cl
  fit$model.matrix <- mm
  fit$response     <- mf[,design$response]
  fit$rand.covar   <- matrix(cbind(rep(1, nrow(mm)), mm[,design$time]), ncol=2)
  fit$ids          <- mf[,design$id]
  fit$n_fixed      <- ncol(mm) - 1

  class(fit) <- "acml"

  if(fit$Code == 3) warning("last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.")
  if(fit$Code == 4) warning("iteration limit exceeded.")
  if(fit$Code == 5) warning("maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.")

  fit
}


#library(MASS, lib.loc="/usr/local/biostat_r/lib/R")
library(MASS)

expit <- function(x) exp(x)/(1+exp(x))
## Standardized Gamma Random effect
gen.std.gam <- function(n.obs, shape, rate){
  (rgamma(n.obs, shape, rate) - shape/rate)/(sqrt(shape/(rate^2)))
}

## Standardized T random effect
gen.std.t <- function(n.obs, df){
  rt(n.obs, df)*sqrt((df-2)/df)
}

## Standardized binormal random effect
gen.std.binormal <- function(n.obs, mn0, mn1, sd0, sd1, p1){
  group <- rbinom(n.obs, 1, p1)
  val <- rnorm(n.obs, mn0, sd0)*(1-group)+rnorm(n.obs, mn1, sd1)*group
  grp.tmp <- rbinom(50000, 1, p1)
  val.tmp <- rnorm(50000, mn0, sd0)*(1-grp.tmp)+rnorm(50000, mn1, sd1)*grp.tmp

  out <- (val-mean(val.tmp))/(sqrt(var(val.tmp)))
  out
}

## Standardized mixture of normals with difference variances
gen.std.mixnorm <- function(n.obs, mn0, mn1, sd0, sd1, p1, group){
  val     <- rnorm(n.obs, mn0, sd0)*(1-group)+rnorm(n.obs, mn1, sd1)*group
  grp.tmp <- rbinom(50000, 1, p1)
  val.tmp <- rnorm(50000, mn0, sd0)*(1-grp.tmp)+rnorm(50000, mn1, sd1)*grp.tmp
  out     <- (val-mean(val.tmp))/(sqrt(var(val.tmp)))
  out
}
GenerateX <- function(N, n, prev.grp, c.parm){
  id   <- rep(1:N, each=n)
  time <- rep(c(0:(n-1)), N)
  grp.tmp <- rbinom(N,1, prev.grp)
  conf.tmp <- rnorm(N, c.parm[1]+grp.tmp*c.parm[2], 1)
  grp  <- rep(grp.tmp, each=n)
  conf <- rep(conf.tmp, each=n)
  out <- data.frame(id=id, time=time, grp=grp, conf=conf)
  out
}

GenerateY <- function(X, Z, id, beta, sig.b0 = 0.25, sig.b1 = 0.25, rho = 0, sig.e = 0.5, RanefDist, ErrorDist){
  lp <- X %*% beta
  cov.mat  <- matrix(c(sig.b0^2, rho*sig.b0*sig.b1, rho*sig.b0*sig.b1, sig.b1^2),2,2)
  ni       <- c(unlist(tapply(id, id, length)))
  N        <- length(unique(id))
  sum.ni   <- length(id)
  if (RanefDist=="Gaussian"){bi <- mvrnorm(N, mu=c(0,0), Sigma=cov.mat)
  }else if (RanefDist=="Gamma5"){b0i <- gen.std.gam(n.obs=N, shape=5, rate=sqrt(3))*sig.b0
  b1i <- gen.std.gam(n.obs=N, shape=5, rate=sqrt(3))*sig.b1
  bi  <- cbind(b0i, b1i)
  }else if (RanefDist=="Binormal1"){b0i <- gen.std.binormal(n.obs=N, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.b0
  b1i <- gen.std.binormal(n.obs=N, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.b1
  bi  <- cbind(b0i, b1i)
  }else if (RanefDist=="MixNorm2"){b0i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.b0
  b1i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.b1
  bi  <- cbind(b0i, b1i)
  }
  b        <- cbind(rep(bi[,1], ni),rep(bi[,2], ni))

  if (ErrorDist=="Gaussian"){error <- mvrnorm(sum.ni, mu=0, Sigma=sig.e)
  }else if (ErrorDist=="Gamma5"){error <- gen.std.gam(n.obs=sum.ni, shape=5, rate=sqrt(3))*sig.e
  }else if (RanefDist=="Binormal1"){error <- gen.std.binormal(n.obs=sum.ni, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.e
  }else if (RanefDist=="MixNorm2"){error <- gen.std.mixnorm(n.obs=sum.ni, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.e
  }

  Y <- X %*% beta + Z[,1]*b[,1] + Z[,2]*b[,2] + error
  return(Y)
}
## generate a missing at random mechanism
GenerateMAR <- function(Y, X, id, param, cutpoint.drop){
  ## set param[2] to 0 for MCAR
  keep <- rep(1, length(Y))
  L    <- length(Y)
  for (j in 2:L){
    if (id[j]==id[j-1] & keep[j-1]==1){ keep[j] <- rbinom(1,1,expit(param[1]+param[2]*I(Y[j-1]<cutpoint.drop) + param[3]*X[j]))
    }else{                              keep[j] <- 1}
  }
  keep
}

## do subject-specific linear regressions. Important to note that data must contain variables Y and time
LinRegFn <- function(data){  X  <- cbind(1, data$time)
Xt <- t(X)
Y  <- data$Y
solve(Xt %*% X) %*% Xt %*% Y}

## calc subject specific intercepts and slopes and output them with the same length as the longitudinal data
CalcSSIntSlp <- function( Y, time, id){
  data.tmp  <- data.frame(id=id, Y=Y, time=time)
  data.list <- split(data.tmp, id)
  L.id      <- c(unlist(tapply(id,id,length)))
  mtx       <- matrix(unlist(lapply(data.list, LinRegFn)), byrow=TRUE, ncol=2)
  out       <- list(Int = rep(mtx[,1], L.id), Slp = rep(mtx[,2], L.id))}

## Calculate the cutpoints for defining sampling strata for the int- slope- and biv-based sampling
## based on the proportion we want in the central region
## If you have a problem with this function it is likely due to the search for the central region
## under bivariate sampling.  You may have to adjust Del.  I think it is due to the discreteness of
## due to insufficient sample size.
est.cutoffs <- function(Y, time, id, PropInCentralRegion){
  p         <- PropInCentralRegion
  data.tmp  <- data.frame(id=id, Y=Y, time=time)
  data.list <- split(data.tmp, id)
  print("Running individual regressions")
  out      <- matrix(unlist(lapply(data.list, LinRegFn)), byrow=TRUE, ncol=2)

  Ints <- quantile(out[,1], c((1-p)/2,(1+p)/2))
  Slps <- quantile(out[,2], c((1-p)/2,(1+p)/2))

  q1 <- .99
  Del <- 1
  while (Del>0.003){
    q1 <- q1-.001
    Del <- abs(mean(out[,1] > quantile(out[,1], probs=1-q1) & out[,1] < quantile(out[,1], probs=q1) &
                      out[,2] > quantile(out[,2], probs=1-q1) & out[,2] < quantile(out[,2], probs=q1)) - PropInCentralRegion)
  }

  out <- list(IntCutUniv = Ints,
              SlpCutUniv = Slps,
              IntCutBiv   = c(quantile(out[,1], probs=1-q1), quantile(out[,1], probs=q1)),
              SlpCutBiv   = c(quantile(out[,2], probs=1-q1), quantile(out[,2], probs=q1)))
  out
}

identify.stratum <- function(Y, time, id, w.function, cutpoints, Int, Slp){
  ## under bivar sampling cutpoints should be of the form (int.lower, int.upper, slp.lower, slp.upper)

  if (w.function == "intercept") stratum <- ifelse( Int <  cutpoints[1], 1,
                                                    ifelse( Int >= cutpoints[2], 3, 2))
  if (w.function == "slope")     stratum <- ifelse( Slp <  cutpoints[1], 1,
                                                    ifelse( Slp >= cutpoints[2], 3, 2))
  if (w.function=="bivar") stratum <- ifelse(Int >= cutpoints[1] & Int < cutpoints[2] & Slp >=cutpoints[3] & Slp <cutpoints[4], 1, 2)
  stratum
}


ods.sampling <- function(id.long,          # id in long format
                         stratum.long,     # stratum identifier in long format)
                         SamplingStrategy, # "IndepODS" or "DepODS"
                         NsPerStratum){    # target number of subjects sampled per stratum.  This is exact if using Dependent sampling
  # This should be of length 3 for univariate sampling and of length 2 for bivariate sampling

  strat.1       <- c(unlist(tapply(stratum.long, id.long, unique)))
  id.1          <- c(unlist(tapply(id.long, id.long, unique)))
  ni            <- c(unlist(tapply(id.long, id.long, length)))
  N             <- length(id.1)
  NPerStratum   <- c(unlist(tapply(id.1, strat.1, length)))
  SampleTooMany <- any(NsPerStratum>NPerStratum)
  if (SampleTooMany){ print("Warning: You want to sample more people than you have in one of the strata.  Sampling from that stratum with probability 1")
    WhichStratum <- which(NsPerStratum>NPerStratum)
    NsPerStratum[WhichStratum] <- NPerStratum[WhichStratum]}
  SampProbs <- NsPerStratum/NPerStratum
  SampProb.1 <- ifelse(strat.1==1, SampProbs[1],
                       ifelse(strat.1==2, SampProbs[2], SampProbs[3]))
  if (SamplingStrategy=="IndepODS"){Samp <- rbinom(N, 1, SampProb.1)}
  if (SamplingStrategy=="DepODS"){  Sampled.ids <- NULL
  for (mm in 1:length(NPerStratum)){
    Sampled.ids <- c( Sampled.ids, c(sample(id.1[strat.1==mm], NsPerStratum[mm], replace=FALSE)))}
  Samp <- ifelse(id.1 %in% Sampled.ids, 1, 0)}
  Sampled <- list(Sampled=rep(Samp, ni),SampProbi=rep(SampProb.1, ni), SampProbs=SampProbs)
  Sampled
}

## note that we really do not need quants, PopnQuants, and w.function but to make the fitting function work
random.sampling <- function(id.long, n=225){
  s <- sample(unique(id.long), n)
  Sampled <- as.integer(id.long %in% s)
  Sampled
}

###############################################################################################
library(lme4)
library(mvtnorm)
library(mitools)


RunMods <- function(params,
                    NsPerStratumUniv,
                    N,
                    prev.grp,
                    n.imp,
                    p.central,
                    conf.param,
                    ni){

  params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5))
  NsPerStratumUniv=c(150,100,150)
  NsPerStratumBiv=c(100,300)

  N=2000
  prev.grp=0.3
  n.imp=50
  p.central=0.8
  conf.param = c(-0.25, 0.5)
  ni = c(4, 6)


  inits            <- params
  NsRand           <- sum(NsPerStratumUniv)
  ##################################################

  ## Generate a full cohort dataframe
  dat      <- GenerateX(N=N, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
  dat$y    <- GenerateY(X=cbind(1, dat$time, dat$grp, dat$conf, dat$time*dat$grp), Z=cbind(1, dat$time), id=dat$id,
                        beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")
  dat      <- dat[order(dat$id, dat$time),]
  droptime <- rep(sample(seq(ni[1], ni[2]), N, replace=TRUE), each=ni[2]) - 1
  dat      <- dat[dat$time<=droptime,]
  cutoffs   <- est.cutoffs(Y=dat$y, time=dat$time, id=dat$id, PropInCentralRegion=p.central)

  ## Calculate subject specific intercepts and slopes
  IntSlps <- CalcSSIntSlp( Y=dat$y, time=dat$time, id=dat$id)
  dat$Int <- IntSlps[[1]]
  dat$Slp <- IntSlps[[2]]

  ## identify stratum membership
  dat$StratInt <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="intercept", cutpoints=cutoffs$IntCutUniv, Int=dat$Int, Slp=dat$Slp)
  dat$StratSlp <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="slope",     cutpoints=cutoffs$SlpCutUniv, Int=dat$Int, Slp=dat$Slp)
  # dat$StratBiv <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="bivar",     cutpoints=c(cutoffs$IntCutBiv, cutoffs$SlpCutBiv), Int=dat$Int, Slp=dat$Slp)

  ## identify the sample along with individual sampling probs and stratum sampling probs
  SampledRan <- random.sampling(id.long=dat$id, n=NsRand)
  SampledInt <- ods.sampling(id.long=dat$id, stratum.long=dat$StratInt, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  SampledSlp <- ods.sampling(id.long=dat$id, stratum.long=dat$StratSlp, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  # SampledBiv <- ods.sampling(id.long=dat$id, stratum.long=dat$StratBiv, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumBiv)

  ## add who will be sampled for each design to dat
  dat$SampledRan <- SampledRan
  dat$SampledInt <- SampledInt[["Sampled"]]
  dat$SampledSlp <- SampledSlp[["Sampled"]]
  # dat$SampledMix1 <- dat$SampledInt*(dat$id <= N/3) + dat$SampledSlp*(dat$id > N/3)
  # dat$SampledMix2 <- dat$SampledInt*(dat$id <= 2*N/3) + dat$SampledSlp*(dat$id > 2*N/3)
  # dat$SampledBiv <-  SampledBiv[["Sampled"]]

  ## Added subject specific sampling weights for each design if doing a weighted likelihood analysis
  dat$WeightsInt <- 1/SampledInt[["SampProbi"]]
  dat$WeightsSlp <- 1/SampledSlp[["SampProbi"]]
  # dat$WeightsMix1 <- 1/(SampledInt[["SampProbi"]]*(dat$id <= N/3) + SampledSlp[["SampProbi"]]*(dat$id > N/3))
  # dat$WeightsMix2 <- 1/(SampledInt[["SampProbi"]]*(dat$id <= 2*N/3) + SampledSlp[["SampProbi"]]*(dat$id > 2*N/3))
  # dat$WeightsBiv <- 1/SampledBiv[["SampProbi"]]
  ## for acml.lmem (not acml.lmem2)
  dat$SampProbiInt <- SampledInt[["SampProbi"]]
  dat$SampProbiSlp <- SampledSlp[["SampProbi"]]
  # dat$SampProbiBiv <- SampledBiv[["SampProbi"]]

  ## This is used if we are not doing sampling weights
  dat$NoWeighting <-1

  dat$IntProbLow  <- SampledInt[["SampProbs"]][1]
  dat$IntProbMid  <- SampledInt[["SampProbs"]][2]
  dat$IntProbHigh <- SampledInt[["SampProbs"]][3]
  dat$SlpProbLow  <- SampledSlp[["SampProbs"]][1]
  dat$SlpProbMid  <- SampledSlp[["SampProbs"]][2]
  dat$SlpProbHigh <- SampledSlp[["SampProbs"]][3]

  dat$IntCutoff1 <- cutoffs$IntCutUniv[1]
  dat$IntCutoff2 <- cutoffs$IntCutUniv[2]
  dat$SlpCutoff1 <- cutoffs$SlpCutUniv[1]
  dat$SlpCutoff2 <- cutoffs$SlpCutUniv[2]

  ## Define w.function within dat for each design
  dat$Int.w <- "intercept"
  dat$Slp.w <- "slope"
  # dat$Mix1.w <- ifelse(dat$id <= N/3, dat$Int.w, dat$Slp.w)
  # dat$Mix2.w <- ifelse(dat$id <= 2*N/3, dat$Int.w, dat$Slp.w)
  # dat$Biv.w <- "bivar"

  ## Stratum Sampling probabilities
  SampProbRan <- c(1,1,1)
  SampProbInt <- SampledInt[["SampProbs"]]
  SampProbSlp <- SampledSlp[["SampProbs"]]
  # SampProbBiv <- SampledBiv[["SampProbs"]]

  ## Datasets for sampled subjects
  datRan  <- dat[dat$SampledRan==1,]
  datInt  <- dat[dat$SampledInt==1,]
  datSlp  <- dat[dat$SampledSlp==1,]
  # datMix1 <- dat[dat$SampledMix1==1,]
  # datMix2 <- dat[dat$SampledMix2==1,]
  # datBiv  <- dat[dat$SampledBiv==1,]

  cutpointsRan=cbind(datRan$IntCutoff1, datRan$IntCutoff2)  ## just need these numbers for the function, not used
  SampProbRan=matrix(1, ncol=3, nrow=length(datRan[,1]))
  w.functionRan=datRan$Int.w                                ## just need these numbers for the function, not used

  cutpointsInt=cbind(datInt$IntCutoff1, datInt$IntCutoff2)
  SampProbInt=cbind(datInt$IntProbLow, datInt$IntProbMid, datInt$IntProbHigh)
  w.functionInt=datInt$Int.w

  cutpointsSlp=cbind(datSlp$SlpCutoff1, datSlp$SlpCutoff2)
  SampProbSlp=cbind(datSlp$SlpProbLow, datSlp$SlpProbMid, datSlp$SlpProbHigh)
  w.functionSlp=datSlp$Slp.w

  design.int <- ods(y~time|id,
                    method = w.functionInt[1],
                    p_sample=SampProbInt[1,],
                    data=dat,
                    quantiles=c(0.1, 0.9))

  #### From our new acml code

  Fit.int <- acml(y~time*grp+conf,
                  design.int,
                  dat[dat$id %in% design.int$sample_ids,],
                  init=inits)



  design.slp <- ods(y~time|id,
                    method = w.functionSlp[1],
                    p_sample=SampProbSlp[1,],
                    data=dat,
                    quantiles=c(1, 0.9))

  #### From our new acml code

  Fit.slp <- acml(y~time*grp+conf,
                  design.slp,
                  dat[dat$id %in% design.slp$sample_ids,],
                  # datSlp,
                  init=inits)

  out <- list(
    # rs.est    = Fit.ran$coefficients,
    # rs.mi.est = Fit.ran.mi[["coefficients"]],

    int.est    = Fit.int$coefficients,
    int.se    = sqrt(diag(Fit.int$covariance)),
    int.rob.se    = sqrt(diag(Fit.int$robcov)),


    slp.est    = Fit.slp$coefficients,
    slp.se    = sqrt(diag(Fit.slp$covariance)),
    slp.rob.se    = sqrt(diag(Fit.slp$robcov))

  )

}

coef_int <- NULL
coef_slp <- NULL
se_int <- NULL
se_slp <- NULL
rob_se_int <- NULL
rob_se_slp <- NULL
for (r in 1:10){
  seed_r = seednum + 10000*r
  set.seed(seed_r)
  N=2000
  params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5))
  fit1_temp = RunMods(params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),-0.25,log(3.5)),
                      NsPerStratumUniv=c(150,100,150),
                      N=2000,
                      prev.grp=0.3,
                      n.imp=50,
                      p.central=0.8,
                      conf.param = c(-0.25, 0.5),
                      ni = c(4, 6))
  coef_int <- rbind(coef_int, fit1_temp$int.est)
  coef_slp <- rbind(coef_slp, fit1_temp$slp.est)
  se_int <- rbind(se_int, fit1_temp$int.se)
  se_slp <- rbind(se_slp, fit1_temp$slp.se)
  rob_se_int <- rbind(rob_se_int, fit1_temp$int.rob.se)
  rob_se_slp <- rbind(rob_se_slp, fit1_temp$slp.rob.se)
}

comb <- data.frame(coef_int) %>%  mutate(type = "coef", design = "int") %>%
  rbind(data.frame(se_int) %>%  mutate(type = "se", design = "int")) %>%
  rbind(data.frame(rob_se_int) %>%  mutate(type = "rob_se", design = "int")) %>%
  rbind(data.frame(coef_slp) %>%  mutate(type = "coef", design = "slp")) %>%
  rbind(data.frame(se_slp) %>%  mutate(type = "se", design = "slp")) %>%
  rbind(data.frame(rob_se_slp) %>%  mutate(type = "rob_se", design = "slp"))


write.csv(comb, paste0("/home/yanb1/Jonathan_projects/acml_sim1/acml_sim1_",seednum, ".csv"))




