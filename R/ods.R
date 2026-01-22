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
#' @seealso \code{\link{ods}}
#'
#' @exportS3Method
#' @importFrom graphics abline hist lines
plot.odsdesign <- function(
    x,
    xlab   = "Intercept",
    ylab   = "Slope",
    main   = format(x$formula),
    sub    = paste(x$method, "design"),
    lwd    = 1,
    lty    = 1,
    cutcol = "black",
    cutlwd = 2,
    cutlty = 1,
    label_strata = TRUE,
    ...
  ) {

    if (!is.null(x$sample_ids)) {
      samp <- as.character(colnames(x$z_i)) %in% x$sample_ids
      col_pts <- ifelse(samp,
                        rgb(0, 0, 0, 0.9),   # sampled
                        rgb(0, 0, 0, 0.15))  # not sampled
    } else {
      col_pts <- if (x$method == "mean") "lightgray" else rgb(0, 0, 0, x$p_sample_i)
    }

    if (x$method == "mean") {
      xlab_hist <- if (xlab == "Intercept") "Mean" else xlab
      hist(x$z_i[1, ], xlab = xlab_hist, main = main, sub = sub,
           col = "lightgray", lwd = lwd, lty = lty, ...)
    } else {
      plot(x$z_i["intercept", ], x$z_i["slope", ],
           xlab = xlab, ylab = ylab,
           main = main, sub = sub,
           col  = col_pts, pch = 3,
           cex = 0.4,
           ...)
    }

    usr <- par("usr")

    if (!(x$method %in% c("bivariate","blup.bivariate"))) {

      cp <- x$cutpoints

      for (i in colnames(cp)) {
        if (i %in% c("mean", "intercept")) {
          abline(v = cp[, i], col = cutcol, lty = cutlty, lwd = cutlwd, ...)
        } else {
          abline(h = cp[, i], col = cutcol, lty = cutlty, lwd = cutlwd, ...)
        }
      }

      if (label_strata) {
        if (any(colnames(cp) %in% c("mean", "intercept"))) {
          name_x <- intersect(c("intercept", "mean"), colnames(cp))[1]
          cx <- sort(cp[, name_x])
          bounds_x <- c(usr[1], cx, usr[2])
          mid_x <- (bounds_x[-1] + bounds_x[-length(bounds_x)]) / 2
          y_top <- usr[4] - 0.05 * (usr[4] - usr[3])

          for (k in seq_along(mid_x)) {
            text(mid_x[k], y_top, labels = bquote(R^.(k)))
          }
        }

        if ("slope" %in% colnames(cp)) {
          cy <- sort(cp[, "slope"])
          bounds_y <- c(usr[3], cy, usr[4])
          mid_y <- (bounds_y[-1] + bounds_y[-length(bounds_y)]) / 2
          x_left <- usr[1] + 0.05 * (usr[2] - usr[1])

          for (k in seq_along(mid_y)) {
            text(x_left, mid_y[k], labels = bquote(R^.(k)))
          }
        }
      }

    } else {

      square <- function(xc, yc) {
        lines(xc[c(1, 1, 2, 2, 1)], yc[c(1, 2, 2, 1, 1)],
              col = cutcol, lty = cutlty, lwd = cutlwd, ...)
      }

      n <- nrow(x$cutpoints)
      for (i in 1:(n / 2)) {
        sel <- c(i, n + 1 - i)
        square(x$cutpoints[sel, 1], x$cutpoints[sel, 2])
      }

      if (label_strata) {
        rx <- range(x$cutpoints[, 1])
        ry <- range(x$cutpoints[, 2])
        text(mean(rx), mean(ry), labels = expression(R^1))
        x_out <- usr[1] + 0.1 * (usr[2] - usr[1])
        y_out <- usr[3] + 0.9 * (usr[4] - usr[3])
        text(x_out, y_out, labels = expression(R^2))
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
    prob_intercept = NULL,
    na.action = getOption('na.action'),
    ProfileCol= NULL   ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
)
{
  # Validate arguments
  coll <- makeAssertCollection()

  ## method could only be 3 types: slope|intercept|bivariate|mean
  assert_character(
    method,
    len           = 1,
    any.missing   = FALSE,
    pattern       = "^(slope|intercept|bivariate|mean)$",
    add           = coll
  )

  ## formula type
  assert_formula(formula, add = coll)

  ## y ~ . structure
  assert_true(
    length(formula) == 3L,
    .var.name = "formula must be two-sided, e.g. Response ~ Month|Patient",
    add = coll
  )

  rhs <- formula[[3L]]

  ## Right hand side must be time | id
  ok_rhs <- is.call(rhs) && identical(rhs[[1L]], as.name("|")) && length(rhs) == 3L

  assert_true(
    ok_rhs,
    .var.name = paste(
      "right-hand side of formula must be of the form time_formula|id,",
      "e.g. Response ~ Month|Patient or Response ~ Month + Age|Patient"
    ),
    add = coll
  )

  reportAssertions(coll)

  lhs       <- formula[[2L]]
  time_term <- rhs[[2L]]
  id_term   <- rhs[[3L]]

  ## time is variable name
  time_is_single <- is.name(time_term)

  if (!time_is_single) {
    # Allow expressions like Month + Age + ...
    tt <- stats::terms(time_term)
    time_labels <- attr(tt, "term.labels")
    assert_true(
      length(time_labels) >= 1L,
      .var.name = "time part must contain at least one variable, e.g. Response ~ Month|Patient or Response ~ Month + Age|Patient",
      add = coll
    )
  } else {
    time_labels <- as.character(time_term)
  }

  ## id is a variable name
  assert_true(
    is.name(id_term),
    .var.name = "id part must be a single variable name, e.g. Response ~ Month|Patient",
    add = coll
  )

  ## method vs time
  ## Response ~ Month|Patient can only be c("intercept", "slope", "bivariate", "mean")
  assert_true(
    method %in% c("intercept", "slope", "bivariate", "mean"),
    .var.name = "when formula is Response ~ Month|Patient, method must be 'intercept', 'slope', 'bivariate', or 'mean'",
    add = coll
  )

  response_name <- as.character(lhs)
  time_name     <- as.character(time_term)
  id_name       <- as.character(id_term)

  ## p_sample must always be present，and it is consistent with method
  if (identical(method, "bivariate")) {
    ## square donut has inner and outer probability
    assert_numeric(
      p_sample,
      lower        = 0,
      upper        = 1,
      len          = 2,
      any.missing  = FALSE,
      add          = coll
    )
  } else {
    ## three regions, low, medium, high
    assert_numeric(
      p_sample,
      lower        = 0,
      upper        = 1,
      len          = 3,
      any.missing  = FALSE,
      add          = coll
    )
  }

  ## only one of quantiles or cutpoints can be specified
  assert_true(
    xor(is.null(quantiles), is.null(cutpoints)),
    .var.name = "only one of quantiles or cutpoints can be specified",
    add = coll
  )

  if (identical(method, "bivariate")) {

    ## bivariate: quantiles = PropInCentralRegion
    if (!is.null(quantiles)) {
      assert_numeric(
        quantiles,
        lower        = 0,
        upper        = 1,
        len          = 1L,
        null.ok      = FALSE,
        any.missing  = FALSE,
        add          = coll
      )
    }

    ## bivariate: cutpoints = IntLow, IntHigh, SlpLow, SlpHigh
    if (!is.null(cutpoints)) {
      assert_numeric(
        cutpoints,
        len          = 4L,
        null.ok      = FALSE,
        any.missing  = FALSE,
        add          = coll
      )
    }

  } else {

    ## univariate: n_c = length(p_sample) - 1
    n_c <- length(p_sample) - 1L

    if (!is.null(quantiles)) {
      assert_numeric(
        quantiles,
        lower        = 0,
        upper        = 1,
        len          = n_c,
        null.ok      = FALSE,
        any.missing  = FALSE,
        add          = coll
      )
    }

    if (!is.null(cutpoints)) {
      assert_numeric(
        cutpoints,
        len          = n_c,
        null.ok      = FALSE,
        any.missing  = FALSE,
        add          = coll
      )
    }
  }

  reportAssertions(coll)

  ## model.frame (similar to lm)
  cl      <- match.call()
  mf_call <- match.call(expand.dots = FALSE)
  m       <- match(c("formula", "data", "subset", "na.action"),
                   names(mf_call), 0L)
  mf_call <- mf_call[c(1L, m)]
  mf_call$drop.unused.levels <- TRUE
  mf_call[[1L]]        <- quote(stats::model.frame)
  mf_call[["formula"]] <- stats::as.formula(gsub("\\|", "+", format(formula)))
  mf <- eval(mf_call, parent.frame())

  ## id
  id_idx <- match(id_name, names(mf))
  if (is.na(id_idx)) {
    stop("id variable '", id_name, "' not found in data/model.frame")
  }
  if (!is.integer(mf[[id_idx]]) && !is.numeric(mf[[id_idx]])) {
    mf[[id_idx]] <- as.numeric(as.factor(mf[[id_idx]]))
  }


  ## construct z_i: mean, intercept, slope

  y_name   <- response_name
  time_idx <- match(time_name, names(mf))
  if (is.na(time_idx)) {
    stop("time variable '", time_name, "' not found in data/model.frame")
  }

  reportAssertions(coll)

  ## Response ~ Month|Patient: regressed on time
  z_fun <- function(x) {
    y  <- x[[y_name]]
    tt <- x[[time_name]]
    fit <- stats::lm(y ~ tt, na.action = na.action)
    c(mean(y), stats::coef(fit))
  }

  z_i <- sapply(split(mf, mf[[id_idx]]), z_fun)
  rownames(z_i) <- c("mean", "intercept", "slope")

  ## z_mf: fixed effect matrix
  z_mf <- stats::model.matrix(
    stats::reformulate(time_name, intercept = TRUE),
    data = mf
  )

  ## create/process cutpoints

  if (is.null(cutpoints)) {

    ## from quantiles and data to calculate cutpoints
    if (identical(method, "bivariate")) {

      ## bivariate: quantiles are PropInCentralRegion
      p0   <- quantiles[1]
      ints <- z_i["intercept", ]
      slps <- z_i["slope",     ]

      q1  <- 0.99
      Del <- 1

      ## marginal joint inside around p0
      while (Del > 0.003 && q1 > 0.5) {
        q1 <- q1 - 0.001

        I_low  <- as.numeric(stats::quantile(ints, probs = 1 - q1))
        I_high <- as.numeric(stats::quantile(ints, probs = q1))
        S_low  <- as.numeric(stats::quantile(slps, probs = 1 - q1))
        S_high <- as.numeric(stats::quantile(slps, probs = q1))

        inside <- (ints > I_low & ints < I_high &
                     slps > S_low & slps < S_high)

        Del <- abs(mean(inside) - p0)
      }

      cutpoints <- rbind(
        low  = c(intercept = I_low,  slope = S_low),
        high = c(intercept = I_high, slope = S_high)
      )

    } else if (identical(method, "mixture")) {

      cutpoints <- apply(z_i, 1, quantile, quantiles)[,c("intercept", "slope"), drop=FALSE]
      cutpoints <- matrix(
        cutpoints,
        nrow=2, byrow=FALSE,
        dimnames=list(1:2, c("intercept", "slope"))
      )

      } else {

      ## univariate: intercept or slope
      n_c     <- length(p_sample) - 1L
      cp_main <- as.numeric(stats::quantile(z_i[method, ], probs = quantiles))
      cutpoints <- matrix(
        cp_main,
        nrow     = n_c,
        ncol     = 1L,
        dimnames = list(seq_len(n_c), method)
      )
    }

  } else {

    ## user defined cutpoints
    if (method %in% c("bivariate","mixture")) {
      ## c(IntLow, IntHigh, SlpLow, SlpHigh)
      cutpoints <- rbind(
        low  = c(intercept = cutpoints[1], slope = cutpoints[3]),
        high = c(intercept = cutpoints[2], slope = cutpoints[4])
      )
    } else {
      n_c <- length(p_sample) - 1L
      cutpoints <- matrix(
        cutpoints,
        nrow     = n_c,
        ncol     = 1L,
        dimnames = list(seq_len(n_c), method)
      )
    }
  }


  ## calculate each subject's sampling probability p_sample_i

  if (identical(method, "bivariate")) {

    ## square donut inside/outside
    ints <- z_i["intercept", ]
    slps <- z_i["slope",     ]

    inside <- (ints >  cutpoints["low",  "intercept"] &
                 ints <  cutpoints["high", "intercept"] &
                 slps >  cutpoints["low",  "slope"] &
                 slps <  cutpoints["high", "slope"])

    p_sample_i <- ifelse(inside, p_sample[2], p_sample[1])
    names(p_sample_i) <- colnames(z_i)

  } else if (identical(method, "mixture")) {

    # for mixture design, the sampling probability should be a vector of 3, which is low, medium, high, the same as univariate.

    p_sample_vec <- data.frame(p_intercept = p_sample[as.numeric(
      cut(z_i['intercept',], c(-Inf, t(cutpoints[,'intercept']), Inf)))],
      p_slope = p_sample[as.numeric(
        cut(z_i['slope',],     c(-Inf, t(cutpoints[,'slope']),     Inf)))])
    rownames(p_sample_vec) <- colnames(z_i)

  } else {

    ## univariate：from cutpoints to have n_c+1 regions
    bounds <- as.numeric(cutpoints[, method])
    p_sample_i <- p_sample[
      as.numeric(
        cut(
          z_i[method, ],
          c(-Inf, bounds, Inf)
        )
      )
    ]
    names(p_sample_i) <- colnames(z_i)
  }

  ## sampling ids
  if (identical(method, "mixture")) {
    # prob_intercept is the proportion of sampled based on intercept.
    sampled_by_intercept = stats::rbinom(nrow(p_sample_vec), 1, prob_intercept) > 0
    p_sample_i = ifelse(sampled_by_intercept, p_sample_vec[,"p_intercept"], p_sample_vec[,"p_slope"])
    smpl <- rownames(p_sample_vec)[((stats::rbinom(nrow(p_sample_vec), 1, p_sample_vec[,"p_intercept"]) > 0)*sampled_by_intercept) | ((stats::rbinom(nrow(p_sample_vec), 1, p_sample_vec[,"p_slope"]) > 0)*(1-sampled_by_intercept))]
  }
  else {
    smpl <- names(p_sample_i)[stats::rbinom(length(p_sample_i), 1, p_sample_i) > 0]
  }

  ## FIXME: adding sample() to sample exact numbers of subjects.

  if (is.null(ProfileCol)) {
    ProfileCol <- NA
  }

  structure(
    list(
      call        = cl,
      formula     = formula,
      model.frame = mf,
      method      = method,
      p_sample    = p_sample,
      p_sample_i  = p_sample_i,
      sample_ids  = smpl,
      response    = response_name,
      time        = time_name,
      id          = id_name,
      quantiles   = quantiles,
      cutpoints   = cutpoints,
      z_i         = z_i,
      z_mf        = z_mf,
      n_rand      = 2,
      ProfileCol  = ProfileCol
    ),
    class = "odsdesign"
  )
}
