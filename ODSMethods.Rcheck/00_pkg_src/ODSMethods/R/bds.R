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

  xx <- matrix(
    rep(NA, 24),
    ncol = 4,
    dimnames = list(
      c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.", "Mean"),
      c(names(object$model.frame)[1], rownames(object$z_i))
    )
  )

  xx[1:5, 1] <- stats::quantile(object$model.frame[, 1], na.rm = TRUE, names = FALSE)
  xx[6, 1]   <- mean(object$model.frame[, 1], na.rm = TRUE)

  xx[1:5, 2] <- stats::quantile(object$z_i["mean", ], na.rm = TRUE, names = FALSE)
  xx[6, 2]   <- mean(object$z_i["mean", ], na.rm = TRUE)

  xx[1:5, 3] <- stats::quantile(object$z_i["intercept", ], na.rm = TRUE, names = FALSE)
  xx[6, 3]   <- mean(object$z_i["intercept", ], na.rm = TRUE)

  xx[1:5, 4] <- stats::quantile(object$z_i["slope", ], na.rm = TRUE, names = FALSE)
  xx[6, 4]   <- mean(object$z_i["slope", ], na.rm = TRUE)

  ans$descriptive <- as.table(xx[c(1, 2, 6, 3:5), ])

  ans$N <- c(nrow(object$model.frame), ncol(object$z_i), sum(object$p_sample_i))
  names(ans$N) <- c("N", names(object$model.frame)[3], "E[N_sample]")

  class(ans) <- "summary.bdsdesign"
  ans
}

#' @exportS3Method
plot.bdsdesign <- function(
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
  plot.odsdesign(
    x = x,
    xlab = xlab,
    ylab = ylab,
    main = main,
    sub = sub,
    lwd = lwd,
    lty = lty,
    cutcol = cutcol,
    cutlwd = cutlwd,
    cutlty = cutlty,
    label_strata = label_strata,
    ...
  )
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
#' odsd <- bds(Response ~ Month|Patient, 'bivariate', p_sample=c(0.25, 1),
#'             data=gbti, quantiles=0.8)
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
    method    = NULL,
    p_sample  = NULL,
    weights   = NULL,
    data      = NULL,
    quantiles = NULL,
    cutpoints = NULL,
    subset    = NULL,
    prob_intercept = NULL,
    method_name  = NULL,
    acml_samp_prob_name = NULL,
    cutpoints_name = NULL,
    sampled_name   = NULL,
    na.action = getOption('na.action'),
    ProfileCol= NULL)  ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
{
  # Validate arguments
  coll <- makeAssertCollection()

  have_main <- !is.null(method) &&
    !is.null(p_sample) &&
    xor(is.null(quantiles), is.null(cutpoints))

  have_alt  <- !is.null(method_name) &&
    !is.null(acml_samp_prob_name) &&
    !is.null(cutpoints_name) &&
    !is.null(sampled_name)

  assert_true(
    have_main || have_alt,
    .var.name = paste(
      "Either supply 'method' and 'p_sample' together with exactly one of",
      "'quantiles' or 'cutpoints', or supply all of",
      "'method_name', 'acml_samp_prob_name', 'cutpoints_name', and 'sampled_name'."
    ),
    add = coll
  )

  ## method could only be 3 types: slope|intercept|bivariate|mean
  if (!is.null(method)){
    assert_character(
      method,
      len           = 1,
      any.missing   = FALSE,
      pattern       = "^(slope|intercept|bivariate|mean|mixture)$",
      add           = coll
    )
    assert_true(
      method %in% c("intercept", "slope", "bivariate", "mean", "mixture"),
      .var.name = "when formula is Response ~ Month|Patient, method must be 'intercept', 'slope', 'bivariate', 'mean' or 'mixture'",
      add = coll
    )
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
  }


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

  # ## id is a variable name
  # assert_true(
  #   is.name(id_term),
  #   .var.name = "id part must be a single variable name, e.g. Response ~ Month|Patient",
  #   add = coll
  # )

  ## method vs time
  ## Response ~ Month|Patient can only be c("intercept", "slope", "bivariate", "mean")


  response_name <- as.character(lhs)
  time_name     <- as.character(time_term)
  id_name       <- as.character(id_term)

  ## p_sample must always be present，and it is consistent with method

  if (!is.null(method)){

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
}
  reportAssertions(coll)

  # ## model.frame (similar to lm)
  cl      <- match.call()
  # mf_call <- match.call(expand.dots = FALSE)
  # m       <- match(c("formula", "data", "subset", "na.action"),
  #                  names(mf_call), 0L)
  # mf_call <- mf_call[c(1L, m)]
  # mf_call$drop.unused.levels <- TRUE
  # mf_call[[1L]]        <- quote(stats::model.frame)
  # mf_call[["formula"]] <- stats::as.formula(gsub("\\|", "+", format(formula)))
  # mf <- eval(mf_call, parent.frame())
  #
  #
  # ## id
  # id_idx <- match(id_name, names(mf))
  # if (is.na(id_idx)) {
  #   stop("id variable '", id_name, "' not found in data/model.frame")
  # }
  # if (!is.integer(mf[[id_idx]]) && !is.numeric(mf[[id_idx]])) {
  #   mf[[id_idx]] <- as.numeric(as.factor(mf[[id_idx]]))
  # }
  #
  #
  # ## construct z_i: mean, intercept, slope
  #
  # y_name   <- response_name
  # time_idx <- match(time_name, names(mf))
  # if (is.na(time_idx)) {
  #   stop("time variable '", time_name, "' not found in data/model.frame")
  # }
  #
  # reportAssertions(coll)
  #
  # ## Response ~ Month|Patient: regressed on time
  # z_fun <- function(x) {
  #   y  <- x[[y_name]]
  #   tt <- x[[time_name]]
  #   fit <- stats::lm(y ~ tt, na.action = na.action)
  #   c(mean(y), stats::coef(fit))
  # }
  #
  # z_i <- sapply(split(mf, mf[[id_idx]]), z_fun)
  # rownames(z_i) <- c("mean", "intercept", "slope")
  # colnames(z_i) <- mf[[id_idx]]
  #
  # ## z_mf: fixed effect matrix
  # z_mf <- stats::model.matrix(
  #   stats::reformulate(time_name, intercept = TRUE),
  #   data = mf
  # )

  ftxt  <- paste(deparse(formula), collapse = "")
  parts <- strsplit(ftxt, "\\|", perl = TRUE)[[1]]
  if (length(parts) != 2) stop("Formula must be like: y ~ time | id")

  left_txt <- trimws(parts[1])    # "y ~ time"
  id_name  <- trimws(parts[2])    # "id"
  f_left   <- stats::as.formula(left_txt)

  y_name    <- all.vars(stats::update(f_left, . ~ 0))[1]
  rhs_names <- all.vars(stats::update(f_left, 0 ~ .))
  if (length(rhs_names) != 1) stop("Left side must have exactly one time covariate, e.g. y ~ Month | id")
  time_name <- rhs_names[1]

  mf_formula <- stats::as.formula(paste(y_name, "~", time_name, "+", id_name))

  mf <- stats::model.frame(
    mf_formula,
    data = data,
    na.action = na.action,
    drop.unused.levels = TRUE
  )

  if (!is.character(mf[[id_name]])) {
    mf[[id_name]] <- as.character(mf[[id_name]])
  }

  ## 3) Construct z_i by subject: mean(y), intercept, slope -----------
  split_list <- split(mf, mf[[id_name]])

  z_fun <- function(df_id) {
    y  <- df_id[[y_name]]
    tt <- df_id[[time_name]]

    fit <- stats::lm(y ~ tt, data = df_id) #, na.action = na.action)
    c(mean = mean(y), intercept = stats::coef(fit)[1], slope = stats::coef(fit)[2])
  }

  z_i <- vapply(split_list, z_fun, FUN.VALUE = c(mean = 0, intercept = 0, slope = 0))
  colnames(z_i) <- names(split_list)

  ## 4) Fixed-effect design matrix for (Intercept + time) -------------
  z_mf <- stats::model.matrix(
    stats::reformulate(time_name, intercept = TRUE),
    data = mf
  )


  ## create/process cutpoints
  smpl <- NULL
  if (!is.null(method)) {
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
          dimnames=list(c("low","high"), c("intercept", "slope"))
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


    ## calculate each subject's sampling probability p_sample_i and define each subject's method, cutpoints, sampling probability

    if (identical(method, "bivariate")) {

      ## square donut inside/outside
      ints <- z_i["intercept", ]
      slps <- z_i["slope",     ]

      inside <- (ints >  cutpoints["low",  "intercept"] &
                   ints <  cutpoints["high", "intercept"] &
                   slps >  cutpoints["low",  "slope"] &
                   slps <  cutpoints["high", "slope"])

      p_sample_i <- ifelse(inside, p_sample[1], p_sample[2])
      names(p_sample_i) <- names(split_list)

      method_i         = rep(method, length(names(split_list)))
      names(method_i)  = names(split_list)

      cutpoints_i      = matrix(rep(c(cutpoints["low",  "intercept"], cutpoints["high",  "intercept"], cutpoints["low",  "slope"], cutpoints["high", "slope"]), length(names(split_list))), ncol = 4, byrow = T)
      rownames(cutpoints_i)  = names(split_list)

      acml_samp_prob_i = matrix(rep(p_sample, length(names(split_list))), ncol = 2, byrow = T)
      rownames(acml_samp_prob_i)  = names(split_list)

    } else if (identical(method, "mixture")) {

      # for mixture design, the sampling probability should be a vector of 3, which is low, medium, high, the same as univariate.

      p_sample_vec <- data.frame(p_intercept = p_sample[as.numeric(
        cut(z_i['intercept',], c(-Inf, t(cutpoints[,'intercept']), Inf)))],
        p_slope = p_sample[as.numeric(
          cut(z_i['slope',],     c(-Inf, t(cutpoints[,'slope']),     Inf)))])
      rownames(p_sample_vec) <- names(split_list)


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
      names(p_sample_i) <- names(split_list)

      method_i         = rep(method, length(names(split_list)))
      names(method_i)  = names(split_list)

      cutpoints_i      = matrix(rep(cutpoints, length(names(split_list))), ncol = 2, byrow = T)
      rownames(cutpoints_i)  = names(split_list)

      acml_samp_prob_i = matrix(rep(p_sample, length(names(split_list))), ncol = 3, byrow = T)
      rownames(acml_samp_prob_i)  = names(split_list)

    }

    ## sampling ids
    if (identical(method, "mixture")) {
      # prob_intercept is the proportion of sampled based on intercept.
      sampled_by_intercept = stats::rbinom(nrow(p_sample_vec), 1, prob_intercept) > 0
      p_sample_i = ifelse(sampled_by_intercept, p_sample_vec[,"p_intercept"], p_sample_vec[,"p_slope"])
      smpl <- rownames(p_sample_vec)[((stats::rbinom(nrow(p_sample_vec), 1, p_sample_vec[,"p_intercept"]) > 0)*sampled_by_intercept) | ((stats::rbinom(nrow(p_sample_vec), 1, p_sample_vec[,"p_slope"]) > 0)*(1-sampled_by_intercept))]

      method_i         = ifelse(sampled_by_intercept, "intercept","slope")
      names(method_i)  = names(split_list)

      cutpoints_i      = matrix(unlist(lapply(sampled_by_intercept, function(i) {c(cutpoints["low",  "intercept"], cutpoints["high",  "intercept"]) * i + (1-i)* c(cutpoints["low",  "slope"], cutpoints["high", "slope"])})),ncol = 2, byrow = T)
      rownames(cutpoints_i)  = names(split_list)

      acml_samp_prob_i = matrix(unlist(lapply(sampled_by_intercept, function(i) {p_sample * i + (1-i)* p_sample})),ncol = 3, byrow = T)  # FIXME: allow 6 p_samples
      rownames(acml_samp_prob_i)  = names(split_list)

    }else {
      smpl <- names(p_sample_i)[stats::rbinom(length(p_sample_i), 1, p_sample_i) > 0]
    }
  }

  if (is.null(ProfileCol)) {
    ProfileCol <- NA
  }

  design_tab                    <- data
  if (is.null(method_name)) {
    design_tab$method_i           <- method_i[as.character(design_tab[, id_name])]
  } else {
    design_tab$method_i         <- data[, method_name]
  }

  if (is.null(acml_samp_prob_name)) {
    design_tab$acml_samp_prob_i   <- lapply(as.character(design_tab[, id_name]), function(i) acml_samp_prob_i[i,])
  } else {
    design_tab$acml_samp_prob_i         <- data[, acml_samp_prob_name]
  }

  if (is.null(cutpoints_name)) {
    design_tab$cutpoints_i   <- lapply(as.character(design_tab[, id_name]), function(i) cutpoints_i[i,])
  } else {
    design_tab$cutpoints_i         <- data[, cutpoints_name]
  }

  if (is.null(sampled_name)) {
    design_tab$sampled         <- as.numeric(as.character(design_tab[, id_name]) %in% smpl)
  } else {
    design_tab$sampled         <- data[, sampled_name]
  }

  if (is.null(smpl)){
    sample_ids    = data[data[,sampled_name] == 1, id_name]
  } else {
    sample_ids    = smpl
  }

  xcol.phase1 <- NULL
  ests.phase1 <- NULL

  # Return design object
  structure(list(
    call        = cl,
    formula     = formula,
    data        = design_tab,
    model.frame = mf,
    weights     = weights,
    method      = method,
    p_sample    = p_sample,
    p_sample_i  = p_sample_i,
    sample_ids  = sample_ids,
    response    = response_name,
    time        = time_name,
    id          = id_name,
    quantiles   = quantiles,
    cutpoints   = cutpoints,
    z_i         = z_i,
    z_mf        = z_mf,
    n_rand      = 2,
    ProfileCol  = ProfileCol,
    xcol.phase1 = xcol.phase1,
    ests.phase1 = ests.phase1
  ),
  class = "bdsdesign"
  )
}


