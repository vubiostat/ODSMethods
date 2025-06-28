  ##############################################################################
 #
# ODSMethods Statistical methods in outcome dependent sampling
#
# Copyright (C) 2025 Shawn P. Garbett, Bailu Yan
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

#' Scatter plot an ODS design
#'
#' Generate a basic scatter plot of an ODS design and overlay the
#' cut points.
#'
#' @param x the ods design object output from \code{\link{ods}}.
#' @param xlab a title for the x axis: see \code{\link{title}}.
#' @param ylab a title for the y axis: see \code{\link{title}}.
#' @param main an overall title for the plot: see \code{\link{title}}.
#' @param sub a subtitle for the plot: see \code{\link{title}}.
#' @param col A specification for the default plotting color. See \code{\link{par}}.
#' @param lwd The line width, a positive number, defaulting to 1. See \code{\link{par}}.
#' @param lty The line type. Line types can either be specified as an integer
#'   (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash,
#'    6=twodash) or as one of the character strings \code{"blank", "solid",
#'    "dashed", "dotted", "dotdash", "longdash", or "twodash"}, where
#'    \code{"blank"} uses ‘invisible lines’ (i.e., does not draw them).
#'    See \code{\link{par}}.
#' @param cutcol A specification for the cut point line plotting color.
#'    Defaults to 'red'. See \code{\link{par}}.
#' @param cutlwd The cut point line width, a positive number, defaulting to 2.
#'    See \code{\link{par}}.
#' @param cutlty The cut point line type. Line types can either be specified as an integer
#'   (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash,
#'    6=twodash) or as one of the character strings \code{"blank", "solid",
#'    "dashed", "dotted", "dotdash", "longdash", or "twodash"}, where
#'    \code{"blank"} uses ‘invisible lines’ (i.e., does not draw them).
#'    See \code{\link{par}}.
#' @param ... Additional arguments past to \code{\link{plot}}, \code{\link{hist}},
#'    or \code{\link{lines}} depending on context.
#' @seealso [ods()]
#'
#' @exportS3Method
#' @importFrom graphics abline hist lines
plot.odsdesign <- function(
  x,
  xlab   = "Intercept",
  ylab   = "Slope",
  main   = format(x$formula),
  sub    = paste(x$method, "design"),
  col    = if(x$method == 'mean') 'lightgray' else 'black',
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
    plot(x$z_i[2,], x$z_i[3,],
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

#' Specify a given design for Outcome Dependent Sampling (ODS)
#'
#' Specify the design of an outcome dependent sampling routine.
#'
#' @param formula (or an object that can be
#'   coerced to that class): a symbolic description of the sampling
#'   model.
#' @param method `character(1)`; A string that specifies 'slope', 'intercept', 'bivariate', or 'mean'.
#' @param p_sample `numeric(n)`; A numeric vector with sampling probabilities
#'   for each region specified by cutpoints.
#' @param data `data.frame`; an optional data frame, list or environment (or
#'   object coercible by as.data.frame to a data frame) containing the variables
#'   in the model. If not found in data, the variables are taken from
#'   environment(formula), typically the environment from which acml is called.
#' @param quantiles `numeric(n-1)`; quantiles to compute cutpoints of interest.
#'   One of quantiles or cutpoints must be specified. Default is NULL.
#' @param cutpoints `numeric(n-1)`; a specification of the cutpoints for sampling. Default is NULL.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process. (See additional details about how this
#'   argument interacts with data-dependent bases in the ‘Details’ section of
#'   the model.frame documentation.)
#' @param weights	an optional vector of weights to be used in the fitting
#'   process. Should be NULL or a numeric vector. If non-NULL, weighted least
#'   squares is used with weights weights (that is, minimizing sum(w*e^2));
#'   otherwise ordinary least squares is used. See also ‘Details’,
#' @param na.action a function which indicates what should happen when the data
#'   contain NAs. The default is set by the na.action setting of options, and is
#'   na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another
#'   possible value is NULL, no action. Value na.exclude can be useful.
#' @param ... additional arguments.
#' @return Returns an ODS design object.
#' @seealso [plot.odsdesign()]
#' @examples
#' data(gbti)
#'
#' odsd <- ods(Response ~ Month|Patient, 'mean', p_sample=c(1, 0.25, 1),
#'             data=gbti, quantiles=c(0.1, 0.9))
#' summary(odsd)
#' plot(odsd)
#'
#' odsd <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
#'             data=gbti, quantiles=c(0.1, 0.9))
#' summary(odsd)
#' plot(odsd)
#'
#' odsd <- ods(Response ~ Month|Patient, 'slope', p_sample=c(1, 0.25, 1),
#'             data=gbti, quantiles=c(0.1, 0.9))
#' summary(odsd)
#' plot(odsd)
#'
#' odsd <- ods(Response ~ Month|Patient, 'bivariate', p_sample=c(1, 0.25, 1),
#'             data=gbti, quantiles=c(0.1, 0.9))
#' summary(odsd)
#' plot(odsd)
#'
#' @export
#'
#' @importFrom checkmate makeAssertCollection
#' @importFrom checkmate assert_formula assert_numeric assert_true assert_character
#' @importFrom checkmate reportAssertions
#' @importFrom stats model.frame terms as.formula quantile coef
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
  ...
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

  # Return design object
  structure(list(
      call        = cl,
      formula     = formula,
      model.frame = mf,
      method      = method,
      p_sample    = p_sample,
      p_sample_i  = p_sample_i,
      response    = names(mf)[1],
      time        = names(mf)[2],
      id          = names(mf)[3],
      quantiles   = quantiles,
      cutpoints   = cutpoints,
      z_i         = z_i,
      n_rand      = 2     # Number of random effects, slope + intercept
    ),
    class="odsdesign"
  )
}
