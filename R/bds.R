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


#' @exportS3Method
#' @importFrom stats quantile
summary.bdsdesign <- function(object, digits = max(3L, getOption("digits")), ...)
{
  ans <- object[c("call", "cutpoints")]
  ans$digits <- digits
  xx <- matrix(rep(NA, 18), ncol=3,dimnames=list(
    c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max." , "Mean"),
    c(names(object$model.frame)[1], rownames(object$z_i))
  ))
  xx[1:5,1] <- quantile(object$model.frame[,1], na.rm=TRUE, names=FALSE)
  xx[6,1]   <- mean(object$model.frame[,1], na.rm=TRUE)
  xx[1:5,2] <- quantile(object$z_i['intercept',], na.rm=TRUE, names=FALSE)
  xx[6,2]   <- mean(object$z_i['intercept',], na.rm=TRUE)
  xx[1:5,3] <- quantile(object$z_i['slope',], na.rm=TRUE, names=FALSE)
  xx[6,3]   <- mean(object$z_i['slope',], na.rm=TRUE)

  ans$descriptive <- as.table(xx[c(1,2,6,3:5),])

  ans$N <- c(nrow(object$model.frame), ncol(object$z_i), sum(object$p_sample_i))
  names(ans$N) <- c("N", names(object$model.frame)[3], "E[N_sample]")

  class(ans) <- "summary.bdsdesign"
  ans
}


#' Specify a given design for BLUP Dependent Sampling (BDS)
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
#' @return Returns an BDS design object.
#' @seealso [plot.bdsdesign()]
#' @examples
#' data(gbti)
#'
#' odsd <- bds(Response ~ Month|Patient, 'mean', p_sample=c(1, 0.25, 1),
#'             data=gbti, quantiles=c(0.1, 0.9))
#' summary(odsd)
#' plot(odsd)
#'
#' odsd <- bds(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
#'             data=gbti, quantiles=c(0.1, 0.9))
#' summary(odsd)
#' plot(odsd)
#'
#' odsd <- bds(Response ~ Month|Patient, 'slope', p_sample=c(1, 0.25, 1),
#'             data=gbti, quantiles=c(0.1, 0.9))
#' summary(odsd)
#' plot(odsd)
#'
#' odsd <- bds(Response ~ Month|Patient, 'bivariate', p_sample=c(1, 0.25, 1),
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
#' @importFrom lme4 lmer
#' @importFrom stats reformulate
bds <- function(
  formula,
  method,
  p_sample,
  data       = NULL,
  quantiles  = NULL,
  cutpoints  = NULL,
  subset     = NULL,
  na.action  = getOption('na.action'),
  ProfileCol = NULL   ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
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

  # at least one random-effects term
  bars <- lme4::findbars(formula)
  assert_true(
    length(bars) >= 1L,
    .var.name = "formula must contain at least one random-effects term, e.g. (Month|Patient)",
    add = coll
  )

  re1 <- bars[[1L]]
  re_lhs <- re1[[2L]]
  re_id  <- re1[[3L]]

  ## id is a variable name
  assert_true(
    is.name(re_id),
    .var.name = "grouping factor in the first random term must be a single variable name, e.g. (Month|Patient)",
    add = coll
  )

  ## response and id name
  response_name <- as.character(formula[[2L]])
  id_name       <- as.character(re_id)

  ## method vs time
  ## Response ~ Month|Patient can only be c("intercept", "slope", "bivariate", "mean")
  assert_true(
    method %in% c("intercept", "slope", "bivariate", "mean"),
    .var.name = "when formula is Response ~ Month|Patient, method must be 'intercept', 'slope', 'bivariate', or 'mean'",
    add = coll
  )

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

  # Duplicate of lm behavior
  cl      <- match.call()
  mf_call <- match.call(expand.dots = FALSE)
  m       <- match(c("formula", "data", "subset", "na.action"),
                   names(mf_call), 0L)
  mf_call <- mf_call[c(1L, m)]
  mf_call$drop.unused.levels <- TRUE
  mf_call[[1L]]        <- quote(stats::model.frame)
  mf_call[["formula"]] <- stats::as.formula(gsub("\\|", "+", format(formula)))
  mf <- eval(mf_call, parent.frame())

  ## Response ~ Month|Patient: regressed on time
  lme_mod = lme4::lmer(formula = formula,
                       data = data,
                       na.action=na.action)
  re_list <- lme4::ranef(lme_mod)
  if (!id_name %in% names(re_list)) {
    id_name_re <- names(re_list)[1L]
  } else {
    id_name_re <- id_name
  }
  re_mat <- re_list[[id_name_re]]   # rows = subjects, cols = random-effects
  re_cols <- colnames(re_mat)

  has_int   <- "(Intercept)" %in% re_cols
  slope_cols <- setdiff(re_cols, "(Intercept)")

  if (method %in% c("slope","bivariate") && length(slope_cols) < 1L) {
    stop("method '", method, "' requires at least one random slope in the first random-effects term")
  }
  if (method %in% c("intercept","bivariate","mean") && !has_int) {
    stop("method '", method, "' requires a random intercept in the first random-effects term")
  }

  time_name <- if (length(slope_cols) >= 1L) slope_cols[1L] else NA_character_

  z_i <- t(re_mat) # from lme4
  rownames(z_i) <- c("intercept", "slope")

  ## z_mf: fixed effect matrix
  z_mf <- stats::model.matrix(
    stats::reformulate(time_name, intercept = TRUE),
    data = mf
  )

  ## get x.phase1

  X <- stats::model.matrix(lme_mod)      # full fixed-effect design matrix
  xcol.phase1 <- colnames(X)

  ## get beta.phase1:
  beta_fixed <- lme4::fixef(lme_mod)

  vc      <- as.matrix(lme4::VarCorr(lme_mod)[[1]])
  sigma_vc  <- sqrt(diag(vc))
  rho     <- vc[1, 2] / (sigma_vc[1] * sigma_vc[2])  #FIXME: to account for more than 1 random/fixed effects
  sigma_e <- sigma(lme_mod)

  ests.phase1 <- c(beta_fixed, log(sigma_vc), log((1+rho)/(1-rho)), log(sigma_e))

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
    if (identical(method, "bivariate")) {
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
  }

  names(p_sample_i) <- colnames(z_i)

  ## sampling ids

  smpl <- names(p_sample_i)[stats::rbinom(length(p_sample_i), 1, p_sample_i) > 0]


  if (is.null(ProfileCol)) {
    ProfileCol <- NA
  }

  # Return design object
  structure(list(
    call        = cl,
    formula     = formula,
    model.frame = mf,
    method      = paste0("blup.",method),
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
    n_rand      = 2,     # Number of random effects, slope + intercept
    ProfileCol  = ProfileCol, ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
    xcol.phase1 = xcol.phase1,  ## only used for blup sampling
    ests.phase1 = ests.phase1  ## only used for blup sampling
    ),
    class=c("bdsdesign","odsdesign")
  )
}



