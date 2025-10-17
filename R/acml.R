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

#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by q design matrix for random effects
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param Keep.liC If FALSE, the function returns the conditional log likelihood across all subjects.  If TRUE, subject specific contributions and exponentiated subject specific ascertainment corrections are returned in a list.
#' @param xcol.phase1 This only applied if doing BLUP-based sampling.  It is the column numbers of the design matrix x that were used in phase 1 to conduct analyses from which BLUP estimates are calculated. e.g. xcol.phase1 = c(1,2,4) if the first second and fourth columns of x were used in phase 1
#' @param ests.phase1 This only applied if doing BLUP-based sampling.  These are the estimates from the phase 1 analysis.  It is assumed that the columns of the design matrix in phase 1 are a subset of those in phase II.  The estimates should be ordered in the following way and appropriately transformed: (beta, log(variance component SDs), FisherZ(correlation parameters in random effects covariance matrix), log(error SDs)).  The transformed variance component SDs and correlations should be ordered the same way they are ordered in the phase II model
#' @return If Keep.liC=FALSE, conditional log likelihood.  If Keep.liC=TRUE, a two-element list that contains subject specific likelihood contributions and exponentiated ascertainment corrections.
#' @export
#'
LogLikeC2 <- function(y, x, z, w.function, id, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, Weights, Keep.liC=FALSE, xcol.phase1, ests.phase1){

  subjectData    <- CreateSubjectData(id=id,y=y,x=x,z=z,Weights=Weights,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function, xcol.phase1=xcol.phase1, ests.phase1=ests.phase1)
  liC.and.logACi <- lapply(subjectData, LogLikeiC2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)

  if (Keep.liC == FALSE){out <- -1*Reduce('+', unlist(liC.and.logACi))[1]  ## sum ss contributions to liC
  }else{ out <- list(liC    = c(unlist(sapply(unlist(liC.and.logACi), function(x) x[1]))), ## ss contributions to liC
                     logACi = c(unlist(sapply(unlist(liC.and.logACi), function(x) x[2]))))} ## ss contributions to ACi
  out
}



#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function.i, cutpoints.i, wi, mui.phase1
#' @param w.function options include "mean" "intercept" "slope" and "bivar"
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
  wi          =  subjectData[["wi"]]
  mui.phase1  =  subjectData[["mui.phase1"]]
  ni          <- length(yi)
  t.zi        <- t(zi)
  ##########
  ##########
  ##########
  ##########
  wi.tmp <- solve(t.zi %*% zi) %*% t.zi
  if (!(w.function %in% c("bivar", "mvints", "mvslps"))){
    if (w.function %in% c("intercept", "intercept1")){ wi <- wi.tmp[1,]
    } else if (w.function %in% c("slope", "slope1")){  wi <- wi.tmp[2,]
    } else if (w.function %in% c("intercept2")){       wi <- wi.tmp[3,]
    } else if (w.function %in% c("slope2")){           wi <- wi.tmp[4,]
    } else if (w.function=="mean"){                    wi <- t(rep(1/ni, ni))
    }
    wi         <- matrix(wi, 1, ni)
    tmp        <- logACi1q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, mui.phase1)
    logACi     <- tmp[["logACi"]]
    liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  }else {
    if (w.function %in% c("bivar")){ wi <- wi.tmp[c(1,2),]
    } else if (w.function %in% c("mvints")){ wi <- wi.tmp[c(1,3),]
    } else if (w.function %in% c("mvslps")){ wi <- wi.tmp[c(2,4),]
    }
    tmp        <- logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, mui.phase1)
    logACi     <- tmp[["logACi"]]
    liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  }
  ##########
  ##########
  ##########
  ##########
  # if (w.function != "bivar"){
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
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function, SampProb, cutpoints, wi, and mui.phase1
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param sigmae std dev of the measurement error distribution
#' @return gradient of the log transformed ascertainment correction under the bivariate sampling design
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
  wi          =  subjectData[["wi"]]
  mui.phase1  =  subjectData[["mui.phase1"]]
  t.zi        <- t(zi)
  #wi          <- solve(t.zi %*% zi) %*% t.zi
  ###########################
  ###########################
  ###########################
  ###########################
  wi.tmp      <- solve(t.zi %*% zi) %*% t.zi
  if (w.function %in% c("bivar")){         wi <- wi.tmp[c(1,2),]
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
    mu_q    = as.vector(wi %*% ( xi %*% new.param[c(beta.index)] - mui.phase1))
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
  mu_q    = as.vector(wi %*% (xi %*% beta- mui.phase1))
  sigma_q <- wi %*% vi %*% t.wi
  ## for some reason the upper and lower triangles to not always equal, so I am forcing
  ## the upper triangle to equal the lower triangle.
  sigma_q <- (sigma_q + t(sigma_q))/2

  (SampProb[1]-SampProb[2])*Deriv / ACi2q(cutpoints, SampProb, mu_q, sigma_q)
}

#' Calculate \text{V_i = Z_i D t(Z_i) + sig_e^2 I_{n_i}}
#' @param zi n_i by q design matrix for the random effects
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return V_i
#' @export
#'
vi.calc = function(zi, sigma.vc, rho.vc, sigma.e){
  SDMat.RE  = diag(sigma.vc)
  ncolzi    = ncol(zi) ## make sure this equals length(sigma.vc)
  nrowzi    = nrow(zi)
  nERRsd    = length(sigma.e)
  b         = matrix(0,ncolzi,ncolzi)
  b[lower.tri(b, diag=FALSE)] = rho.vc
  CorMat.RE = t(b)+b+diag(rep(1,ncolzi))
  CovMat.RE = SDMat.RE %*% CorMat.RE %*% SDMat.RE
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
ACi1q = function(cutpoints, SampProb, mu_q, sigma_q){
  CDFs = pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
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
ACi2q = function(cutpoints, SampProb, mu_q, sigma_q){
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
#' @param mui.phase1 linear predictor from phase 1 that is used to calculate mu_q if we are doing BLUP based sampling.  If doing ODS then this is a vector of 0s
#' @return log transformed ascertainment correction
#' @export
logACi1q = function(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, mui.phase1){
  vi      = vi.calc(zi, sigma.vc, rho.vc, sigma.e)
  mu      = xi %*% beta
  mu_q    = (wi %*% (mu-mui.phase1))[,1]
  sigma_q = sqrt((wi %*% vi %*% t(wi))[1,1])
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
#' @param mui.phase1 linear predictor from phase 1 that is used to calculate mu_q if we are doing BLUP based sampling.  If doing ODS then this is a vector of 0s
#' @return log transformed ascertainment correction
#' @export
logACi2q = function(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, mui.phase1){
  vi      = vi.calc(zi, sigma.vc, rho.vc, sigma.e)
  mu      = xi %*% beta
  mu_q    = as.vector(wi %*% (mu - mui.phase1))
  sigma_q = wi %*% vi %*% t(wi)
  ## for some reason the upper and lower triangles to not always equal, so I am taking their
  ## average.  Not sure this is a problem here
  ## but doing this to be safe.  Maybe can remove later once understood.
  sigma_q = (sigma_q + t(sigma_q))/2
  #sigma_q[upper.tri(sigma_q)]  = t(sigma_q)[upper.tri(sigma_q)]
  #sigma_q[2,1] = sigma_q[1,2]
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
li.lme = function(yi, xi, beta, vi){
  resid = yi - xi %*% beta
  -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}

#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by q design matrix for random effects
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param Keep.liC If FALSE, the function returns the conditional log likelihood across all subjects.  If TRUE, subject specific contributions and exponentiated subject specific ascertainment corrections are returned in a list.
#' @param xcol.phase1 This only applied if doing BLUP-based sampling.  It is the column numbers of the design matrix x that were used in phase 1 to conduct analyses from which BLUP estimates are calculated. e.g. xcol.phase1 = c(1,2,4) if the first second and fourth columns of x were used in phase 1
#' @param ests.phase1 This only applied if doing BLUP-based sampling.  These are the estimates from the phase 1 analysis.  It is assumed that the columns of the design matrix in phase 1 are a subset of those in phase II.  The estimates should be ordered in the following way and appropriately transformed: (beta, log(variance component SDs), FisherZ(correlation parameters in random effects covariance matrix), log(error SDs)).  The transformed variance component SDs and correlations should be ordered the same way they are ordered in the phase II model
#' @return If Keep.liC=FALSE, conditional log likelihood.  If Keep.liC=TRUE, a two-element list that contains subject specific likelihood contributions and exponentiated ascertainment corrections.
#' @export
#'
LogLikeC2 = function(y, x, z, w.function, id, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, Weights, Keep.liC=FALSE, xcol.phase1, ests.phase1){

  subjectData    = CreateSubjectData(id=id,y=y,x=x,z=z,Weights=Weights,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function, xcol.phase1=xcol.phase1, ests.phase1=ests.phase1)
  liC.and.logACi = lapply(subjectData, LogLikeiC2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)

  if (Keep.liC == FALSE){out = -1*Reduce('+', unlist(liC.and.logACi))[1]  ## sum ss contributions to liC
  }else{ out = list(liC    = c(unlist(sapply(unlist(liC.and.logACi), function(x) x[1]))), ## ss contributions to liC
                    logACi = c(unlist(sapply(unlist(liC.and.logACi), function(x) x[2]))))} ## ss contributions to ACi
  out
}



#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function.i, cutpoints.i, wi, mui.phase1
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
LogLikeiC2 = function(subjectData, beta, sigma.vc, rho.vc, sigma.e){
  yi          = subjectData[["yi"]]
  xi          = subjectData[["xi"]]
  zi          = subjectData[["zi"]]
  Weights.i   = subjectData[["Weights.i"]]
  w.function  = subjectData[["w.function.i"]]
  SampProb    = subjectData[["SampProb.i"]]
  cutpoints   = subjectData[["cutpoints.i"]]
  wi          =  subjectData[["wi"]]
  mui.phase1  =  subjectData[["mui.phase1"]]

  ni          = length(yi)
  t.zi        = t(zi)
  ##########
  ##########
  ##########
  ##########
  #wi.tmp = solve(t.zi %*% zi) %*% t.zi
  if (!(w.function %in% c("bivar", "mvints", "mvslps"))){
    # if (w.function %in% c("intercept", "intercept1")){ wi = wi.tmp[1,]
    # } else if (w.function %in% c("slope", "slope1")){  wi = wi.tmp[2,]
    # } else if (w.function %in% c("intercept2")){       wi = wi.tmp[3,]
    # } else if (w.function %in% c("slope2")){           wi = wi.tmp[4,]
    # } else if (w.function=="mean"){                    wi = t(rep(1/ni, ni))
    # }
    wi         = matrix(wi, 1, ni)
    tmp        = logACi1q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, mui.phase1)
    logACi     = tmp[["logACi"]]
    liC        = li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  }else {
    # if (w.function %in% c("bivar")){ wi = wi.tmp[c(1,2),]
    # } else if (w.function %in% c("mvints")){ wi = wi.tmp[c(1,3),]
    # } else if (w.function %in% c("mvslps")){ wi = wi.tmp[c(2,4),]
    # }
    tmp        = logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, mui.phase1)
    logACi     = tmp[["logACi"]]
    liC        = li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  }
  ##########
  ##########
  ##########
  ##########
  # if (w.function != "bivar"){
  #     if (w.function %in% c("intercept", "intercept1")){ wi= (solve(t.zi %*% zi) %*% t.zi)[1,]
  #     } else if (w.function %in% c("slope", "slope1")){     wi= (solve(t.zi %*% zi) %*% t.zi)[2,]
  #     } else if (w.function %in% c("intercept2")){ wi= (solve(t.zi %*% zi) %*% t.zi)[3,]
  #     } else if (w.function %in% c("slope2")){     wi= (solve(t.zi %*% zi) %*% t.zi)[4,]
  #     } else if (w.function=="mean"){     wi = t(rep(1/ni, ni))
  #     }
  #     wi         = matrix(wi, 1, ni)
  #     tmp        = logACi1q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
  #     logACi     = tmp[["logACi"]]
  #     liC        = li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  #
  # }else{
  #     wi         = solve(t.zi %*% zi) %*% t.zi
  #     tmp        = logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
  #     logACi     = tmp[["logACi"]]
  #     liC        = li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
  #
  # }
  return(c(liC, logACi))
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
  wi          =  subjectData[["wi"]]
  mui.phase1  =  subjectData[["mui.phase1"]]
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
  mu_q      = (wi %*% (mu-mui.phase1))[,1]
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


#' Calculate the gradient of the conditional likelihood for the univariate and bivariate sampling cases across all subjects (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if CheeseCalc=TRUE.
#'
#' Calculate the gradient of the conditional likelihood for the univariate and bivariate sampling cases across all subjects (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if CheeseCalc=TRUE.
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[1,])}, "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[1,])}. "intercept2" (second intercept of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[3,])}, "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[4,])} "bivar" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[c(1,2),])} "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])} "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])}.  NOTE: We also have the same designs but for BLUP based sampling, in which case, the character string should begin with "blup.".  For example "blup.intercept". There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param CheeseCalc If FALSE, the function returns the gradient of the conditional log likelihood across all subjects.  If TRUE, the cheese part of the sandwich esitmator is calculated.
#' @param xcol.phase1 This only applied if doing BLUP-based sampling.  It is the column numbers of the design matrix x that were used in phase 1 to conduct analyses from which BLUP estimates are calculated. e.g. xcol.phase1 = c(1,2,4) if the first second and fourth columns of x were used in phase 1
#' @param ests.phase1 This only applied if doing BLUP-based sampling.  These are the estimates from the phase 1 analysis.  It is assumed that the columns of the design matrix in phase 1 are a subset of those in phase II.  The estimates should be ordered in the following way and appropriately transformed: (beta, log(variance component SDs), FisherZ(correlation parameters in random effects covariance matrix), log(error SDs)).  The transformed variance component SDs and correlations should be ordered the same way they are ordered in the phase II model
#' @return If CheeseCalc=FALSE, gradient of conditional log likelihood.  If CheeseCalc=TRUE, the cheese part of the sandwich estimator is calculated.
#' @export
LogLikeC.Score2 <- function(y, x, z, w.function, id, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, Weights, CheeseCalc=FALSE, xcol.phase1, ests.phase1){
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

  subjectData <- CreateSubjectData(id=id,y=y,x=x,z=z,Weights=Weights,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function, xcol.phase1=xcol.phase1, ests.phase1=ests.phase)

  UncorrectedScorei <- lapply(subjectData, li.lme.score2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)
  Gradienti         <- lapply(UncorrectedScorei, function(x) x[['gr']]) ## create a list of ss contributions to gradient
  UncorrectedScore  <- Reduce('+', Gradienti)  ## Note if using IPW this is actually a corrected score (corrected by the IPW)
  #print(c("blah1", UncorrectedScore))
  ## NOTE HERE: I used the first element of w.function in this call.  This means, for now, we cannot mix bivar with other
  ## sampling schemes.  This also applies to the cheese calculation
  if (!(w.function[[1]] %in% c("bivar","mvints","mvslps"))){
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
    if (!(w.function[[1]] %in% c("bivar","mvints","mvslps"))){ GradiMat <- mapply("+", Gradienti, logACi.Score)
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
#'
#' @param params parameter vector c(beta, log(sigma0), log(sigma1), rho, sigmae)
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param id sum(n_i) vector of subject ids
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[1,])}, "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[1,])}. "intercept2" (second intercept of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[3,])}, "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[2,])}, "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[4,])} "bivar" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects \text{(solve(t.zi %*% zi) %*% t.zi)[c(1,2),])} "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])} "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) \text{solve(t.zi %*% zi) %*% t.zi)[c(1,3),])}.  NOTE: We also have the same designs but for BLUP based sampling, in which case, the character string should begin with "blup.".  For example "blup.intercept". There should be one unique value per subject
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions (bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high).  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maimizing the log likelihood for all other parameters
#'                   while fixing these columns at the values of params[ProfileCol]
#' @param xcol.phase1 This only applied if doing BLUP-based sampling.  It is the column numbers of the design matrix x that were used in phase 1 to conduct analyses from which BLUP estimates are calculated. e.g. xcol.phase1 = c(1,2,4) if the first second and fourth columns of x were used in phase 1
#' @param ests.phase1 This only applied if doing BLUP-based sampling.  These are the estimates from the phase 1 analysis.  It is assumed that the columns of the design matrix in phase 1 are a subset of those in phase II.  The estimates should be ordered in the following way and appropriately transformed: (beta, log(variance component SDs), FisherZ(correlation parameters in random effects covariance matrix), log(error SDs)).  The transformed variance component SDs and correlations should be ordered the same way they are ordered in the phase II model
#' @param Keep.liC If TRUE outputs subject specific conditional log lileihoods to be used for the imputation procedure described in the AOAS paper keep z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @return The conditional log likelihood with a "gradient" attribute (if Keep.liC=FALSE) and subject specific contributions to the conditional likelihood if Keep.liC=TRUE).
#' @export
LogLikeCAndScore2 <- function(params, y, x, z, id, w.function, cutpoints, SampProb, Weights, ProfileCol=NA, Keep.liC=FALSE, xcol.phase1, ests.phase1){
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

  out     <- LogLikeC2( y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e, cutpoints=cutpoints, SampProb=SampProb, Weights=Weights, Keep.liC=Keep.liC, xcol.phase1=xcol.phase1, ests.phase1=ests.phase1)
  GRAD    <- LogLikeC.Score2(y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e, cutpoints=cutpoints, SampProb=SampProb, Weights=Weights, xcol.phase1=xcol.phase1, ests.phase1=ests.phase1)
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
##zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)[3,]), "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)[2,]), "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)[2,]), "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)[4,]) "bivar" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)[c(1,2),]) "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) * t.zi)[c(1,3),]) "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)[c(1,3),]).  There should be one unique value per subject.  NOTE: We also have the same designs but for BLUP based sampling, in which case, the character string should begin with "blup.".  For example "blup.intercept". There should be one unique value per subject but n_i replicates of that value.  Note that w.function should NOT be in the dat dataframe.
#' @param InitVals starting values for c(beta, log(sigma0), log(sigma1), log((1+rho)/(1-rho)), log(sigmae))
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions for individual subjects.  If using a low, medium, high, sampling scheme, this is a sum(n_i) by 2 matrix that must be a distinct object not contained in the dat dataframe.  Each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high.  If using a square doughnut design this should be sum(n_i) by 4 matrix (var1lower, var1upper, var2lower, var2upper). Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region. For low medium high sampling, each row is a vector of length 3 with sampling probabilities for each region. For bivariate stratum sampling each row is a vector of length 2 with sampling probabilities for the inner and outer strata. Each subject should have n_i rows of the same values.  Not in data.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  This should be a variable in the data dataframe. It should only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar.  If not doing IPWL, this must be a vectors of 1s.
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maximizing the log likelihood for all other parameters
#'                   while fixing these columns at the values of InitVals[ProfileCol]
#' @param xcol.phase1 This only applied if doing BLUP-based sampling.  It is the column numbers of the design matrix x that were used in phase 1 to conduct analyses from which BLUP estimates are calculated. e.g. xcol.phase1 = c(1,2,4) if the first second and fourth columns of x were used in phase 1
#' @param ests.phase1 This only applied if doing BLUP-based sampling.  These are the estimates from the phase 1 analysis.  It is assumed that the columns of the design matrix in phase 1 are a subset of those in phase II.  The estimates should be ordered in the following way and appropriately transformed: (beta, log(variance component SDs), FisherZ(correlation parameters in random effects covariance matrix), log(error SDs)).  The transformed variance component SDs and correlations should be ordered the same way they are ordered in the phase II model
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

  id0   = design$id
  id    = data$id = data[ , id0 ]

  fixed.f = model.frame(formula, data)
  fixed.t = attr(fixed.f, "terms")
  y       = model.response(fixed.f,'numeric')
  uy      = unique(y)
  x       = model.matrix(formula, fixed.f)
  z_id    = data.frame(id = design$model.frame[,design$id], design$z_mf)
  z       = as.matrix(merge(data.frame(id = unique(id)), z_id, by = "id", all.x = TRUE))[,-1]

  #if (is.na(SampProb[1])) SampProb = c(1,1,1)
  if (is.null(design$weights)){
    Weights    = data$Weights = rep(1, nrow(data))    #for non-IPW case
  }else{
    Weights0   = design$weights
    Weights    = data$Weights = data[ , Weights0]
  }

  cutpoints   <- rep(list(asplit(design$cutpoints,2)[[1]]), length(y))
  SampProb    <- rep(list(design$p_sample), length(y))
  w.function  <- rep(list(design$method), length(y))

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
                  xcol.phase1=design$xcol.phase1,
                  ests.phase1=design$ests.phase1,
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
                                         ProfileCol=design$ProfileCol,
                                         xcol.phase1=design$xcol.phase1,
                                         ests.phase1=design$ests.phase1)
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
                            CheeseCalc=TRUE,
                            xcol.phase1=design$xcol.phase1,
                            ests.phase1=design$ests.phase1)

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

#' Create a list of subject-specific data
#'
#' @param id sum(n_i) vector of subject ids
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param xcol.phase1 This only applied if doing BLUP-based sampling.  It is the column numbers of the design matrix x that were used in phase 1 to conduct analyses from which BLUP estimates are calculated. e.g. xcol.phase1 = c(1,2,4) if the first second and fourth columns of x were used in phase 1
#' @param ests.phase1 This only applied if doing BLUP-based sampling.  These are the estimates from the phase 1 analysis.  It is assumed that the columns of the design matrix in phase 1 are a subset of those in phase II.  The estimates should be ordered in the following way and appropriately transformed: (beta, log(variance component SDs), FisherZ(correlation parameters in random effects covariance matrix), log(error SDs)).  The transformed variance component SDs and correlations should be ordered the same way they are ordered in the phase II model
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @export
CreateSubjectData <- function(id,y,x,z,Weights,SampProb,cutpoints,w.function, xcol.phase1, ests.phase1){

  if (!is.null(xcol.phase1)){
    npar.phase1   = length(ests.phase1)
    x.phase1      = x[,xcol.phase1]
    nbeta.phase1  = ncol(x.phase1)
    nVCsd.phase1  = ncol(z)
    nVCrho.phase1 = choose(nVCsd.phase1,2)
    nERRsd.phase1 = npar.phase1-nbeta.phase1-nVCsd.phase1-nVCrho.phase1

    beta.index.phase1   = c(1:nbeta.phase1)
    vc.sd.index.phase1  = nbeta.phase1 + (c(1:nVCsd.phase1))
    vc.rho.index.phase1 = nbeta.phase1 + nVCsd.phase1 + (c(1:nVCrho.phase1))
    err.sd.index.phase1 = nbeta.phase1 + nVCsd.phase1 + nVCrho.phase1 + c(1:nERRsd.phase1)

    beta.phase1     = ests.phase1[beta.index.phase1]
    sigma.vc.phase1 = exp(ests.phase1[vc.sd.index.phase1])
    rho.vc.phase1   = (exp(ests.phase1[vc.rho.index.phase1])-1) / (exp(ests.phase1[vc.rho.index.phase1])+1)
    sigma.e.phase1  = exp(ests.phase1[err.sd.index.phase1])

    D.tmp         = diag(sigma.vc.phase1)
    b             = matrix(0,nVCsd.phase1,nVCsd.phase1)
    b[lower.tri(b, diag=FALSE)] = rho.vc.phase1
    R.tmp         = t(b)+b+diag(rep(1,nVCsd.phase1))
    D.phase1      = D.tmp %*% R.tmp %*% D.tmp
    x.phase1.tmp  = split(x.phase1, id)
    ncol.x.phase1 = ncol(x.phase1)
  }

  id.tmp          <- split(id,id)
  y.tmp           <- split(y,id)
  x.tmp           <- split(x,id)
  z.tmp           <- split(z,id)
  Weights.tmp     <- split(Weights,id)
  SampProb.tmp    <- split(SampProb,id)
  cutpoints.tmp   <- split(cutpoints,id)
  w.function.tmp  <- split(w.function,id)

  ncol.x         = ncol(x)
  ncol.z         = ncol(z)
  ncol.SampProb  = ncol(SampProb)
  ncol.cutpoints = ncol(cutpoints)

  subjectData <- vector('list', length=length(unique(id)))
  subjectData <- list()
  uid <- as.character(unique(id))
  for(j in seq(along=uid)){
    i          <- uid[j]
    zi          = matrix(z.tmp[[i]], ncol=ncol.z)
    xi          = matrix(x.tmp[[i]], ncol=ncol.x)
    yi          = y.tmp[[i]]
    w.functioni = unique(w.function.tmp[[i]])
    t.zi        = t(zi)
    ni          = length(yi)
    ######### Calculate wi and mui.phase1 under blup sampling schemes
    if (!is.null(xcol.phase1)){
      xi.phase1  = matrix(x.phase1.tmp[[i]], ncol=ncol.x.phase1)
      Vi.phase1  = vi.calc(zi=zi, sigma.vc=sigma.vc.phase1, rho.vc=rho.vc.phase1, sigma.e=sigma.e.phase1)
      wi.tmp     = D.phase1 %*% t.zi %*% solve(Vi.phase1)
      mui.phase1 = xi.phase1 %*% beta.phase1
      ####################################################
      if (!(w.functioni %in% c("blup.bivar", "blup.mvints", "blup.mvslps"))){
        if (w.functioni %in% c("blup.intercept", "blup.intercept1")){ wi = wi.tmp[1,]
        } else if (w.functioni %in% c("blup.slope", "blup.slope1")){  wi = wi.tmp[2,]
        } else if (w.functioni %in% c("blup.intercept2")){            wi = wi.tmp[3,]
        } else if (w.functioni %in% c("blup.slope2")){                wi = wi.tmp[4,]
        }
        wi         = matrix(wi, 1, ni)
      }else if (w.functioni %in% c("blup.bivar")){ wi = wi.tmp[c(1,2),]
      } else if (w.functioni %in% c("blup.mvints")){ wi = wi.tmp[c(1,3),]
      } else if (w.functioni %in% c("blup.mvslps")){ wi = wi.tmp[c(2,4),]
      } else {stop("You have not chosen an appropriate w.function for BLUP sampling")
      }
      ####################################################
      ######### Calculate wi and mui.phase1 (which equals rep(0,nrow(zi))) under ODS schemes
    }else{
      wi.tmp     = solve(t.zi %*% zi) %*% t.zi
      mui.phase1 = rep(0, nrow(xi))
      ####################################################
      if (!(w.functioni %in% c("bivariate", "mvints", "mvslps"))){
        if (w.functioni %in% c("intercept", "intercept1")){ wi = wi.tmp[1,]
        } else if (w.functioni %in% c("slope", "slope1")){  wi = wi.tmp[2,]
        } else if (w.functioni %in% c("intercept2")){       wi = wi.tmp[3,]
        } else if (w.functioni %in% c("slope2")){           wi = wi.tmp[4,]
        } else if (w.functioni=="mean"){                    wi = t(rep(1/ni, ni))
        }
        wi         = matrix(wi, 1, ni)
      }else if (w.functioni %in% c("bivariate")){ wi = wi.tmp[c(1,2),]
      } else if (w.functioni %in% c("mvints")){ wi = wi.tmp[c(1,3),]
      } else if (w.functioni %in% c("mvslps")){ wi = wi.tmp[c(2,4),]
      } else {stop("You have not chosen an appropriate w.function for ODS sampling")
      }
      ####################################################
    }
    subjectData[[j]] <- list(idi          = as.character(unique(id.tmp[[i]])),
                             xi           = xi,
                             zi           = zi,
                             yi           = yi,
                             wi           = wi,
                             mui.phase1   = mui.phase1,
                             Weights.i    = unique(Weights.tmp[[i]]),
                             SampProb.i   = matrix(SampProb.tmp[[i]][[1]], ncol=length(SampProb.tmp[[i]][[1]])), # assuming only one sampprob for each subject.
                             w.function.i = as.character(w.functioni),
                             cutpoints.i  = matrix(cutpoints.tmp[[i]][[1]], ncol=length(cutpoints.tmp[[i]][[1]])))
  }
  names(subjectData) <- uid
  subjectData
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

  mm <- model.matrix(formula, data)

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

