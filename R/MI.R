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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>

#' @importFrom stats cov
#' @noRd
.mi_rubin_multivariate <- function(coef_mi, cov_mi) {
  M <- nrow(coef_mi)
  beta_bar <- colMeans(coef_mi)
  U_bar <- Reduce("+", cov_mi) / M
  B <- cov(coef_mi)
  T_mi <- U_bar + (1 + 1 / M) * B
  dimnames(T_mi) <- list(colnames(coef_mi), colnames(coef_mi))
  dimnames(U_bar) <- list(colnames(coef_mi), colnames(coef_mi))
  dimnames(B) <- list(colnames(coef_mi), colnames(coef_mi))

  list(beta = beta_bar, vcov = T_mi, U_bar = U_bar, B = B)
}

.mi_design_cutpoints_vector <- function(design) {
  if (is.null(design$cutpoints)) stop("The design object must contain cutpoints.")

  if (identical(design$method, "bivariate")) {
    return(c(
      design$cutpoints["low", "intercept"],
      design$cutpoints["high", "intercept"],
      design$cutpoints["low", "slope"],
      design$cutpoints["high", "slope"]
    ))
  }

  as.numeric(design$cutpoints)
}

.mi_build_design_like <- function(design, data, weights_name, sampled_name) {
  if (is.null(design$method) || is.null(design$p_sample)) {
    stop("MI currently requires an ods design created with method and p_sample.")
  }
  if (identical(design$method, "mixture")) {
    stop("MI currently supports mean, intercept, slope, and bivariate ODS designs.")
  }

  ods(
    formula = design$formula,
    method = design$method,
    p_sample = rep(1, length(design$p_sample)),
    data = data,
    cutpoints = .mi_design_cutpoints_vector(design),
    weights = weights_name,
    sampled_name = sampled_name,
    ProfileCol = design$ProfileCol
  )
}

.mi_formula_rhs <- function(formula) {
  if (!inherits(formula, "formula")) stop("exposure_model must be a formula.")
  if (length(formula) == 2L) return(paste(deparse(formula[[2L]]), collapse = " "))
  if (length(formula) == 3L) return(paste(deparse(formula[[3L]]), collapse = " "))
  stop("exposure_model must be one-sided or two-sided formula.")
}

.mi_default_exposure_model <- function(formula, exposure_name, time_name) {
  vars <- all.vars(formula[[3L]])
  vars <- setdiff(vars, c(exposure_name, time_name))
  if (length(vars) == 0) return(stats::as.formula("~ 1"))
  stats::reformulate(vars)
}

.mi_coef_names <- function(fit) {
  fixed_names <- colnames(fit$model.matrix)
  n_rand <- fit$design$n_rand
  if (n_rand == 2L) {
    random_names <- c(
      "log_sd_(Intercept)",
      paste0("log_sd_", fit$design$time),
      "2_atanh_rho",
      "log_sd_residual"
    )
  } else {
    random_names <- c(
      paste0("log_sd_re", seq_len(n_rand)),
      paste0("2_atanh_rho", seq_len(choose(n_rand, 2L))),
      "log_sd_residual"
    )
  }
  c(fixed_names, random_names)
}

# Type-1 PMM working model: donors use alpha_hat, recipients use alpha_draw.
.mi_pmm_model <- function(H_obs, X_obs, H_mis, ridge)
{
  S <- crossprod(H_obs)
  if (any(diag(S) == 0)) stop("The PMM model contains an all-zero predictor.")
  V <- solve(S + diag(diag(S) * ridge, nrow = ncol(S)))
  alpha_hat <- as.numeric(V %*% crossprod(H_obs, X_obs))
  names(alpha_hat) <- colnames(H_obs)
  pred_donor <- as.numeric(H_obs %*% alpha_hat)
  RSS <- sum((X_obs - pred_donor)^2)
  df_resid <- nrow(H_obs) - ncol(H_obs)
  if (df_resid <= 0) stop("The PMM model needs more donors than coefficients.")
  if (!is.finite(RSS) || RSS <= 0) stop("Invalid PMM RSS or residual degrees of freedom.")

  list(
    coefficients = alpha_hat,
    V            = V,
    V_half       = chol(V),
    pred_donor   = pred_donor,
    RSS          = RSS,
    df_resid     = df_resid,
    H_mis        = H_mis,
    X_obs        = X_obs
  )
}

.mi_pmm_draw <- function(exposure_fit, n_donors)
{
  sigma2_draw <- exposure_fit$RSS / stats::rchisq(1, df = exposure_fit$df_resid)
  z1 <- stats::rnorm(length(exposure_fit$coefficients))
  alpha_draw <- exposure_fit$coefficients +
    sqrt(sigma2_draw) * as.numeric(z1 %*% exposure_fit$V_half)
  pred_recipient <- as.numeric(exposure_fit$H_mis %*% alpha_draw)
  d_use <- min(n_donors, length(exposure_fit$X_obs))
  exposure_imp <- exposure_fit$X_obs[rep(NA_integer_, length(pred_recipient))]

  for (i in seq_along(pred_recipient)) {
    distance_i <- abs(exposure_fit$pred_donor - pred_recipient[i])
    nearest_i <- order(distance_i)[seq_len(d_use)]
    selected_position <- sample.int(length(nearest_i), size = 1L)
    exposure_imp[i] <- exposure_fit$X_obs[nearest_i[selected_position]]
  }
  exposure_imp
}

MI_internal <- function(formula, design, exposure, exposure_model, M, n_donors, ridge, data, sampled_col, na.action, verbose, init, keep_completed)
{
  id_name <- design$id
  time_name <- design$time
  response_name <- design$response
  needed <- c(id_name, time_name, response_name, exposure, all.vars(formula))
  if (length(setdiff(needed, names(data))) > 0L) {
    stop("Missing required column(s): ", paste(setdiff(needed, names(data)), collapse = ", "))
  }
  if (!identical(all.vars(formula[[2L]]), response_name)) {
    stop("must have same response variable as design")
  }
  if (!(is.numeric(data[[exposure]]) || is.logical(data[[exposure]]) ||
        is.factor(data[[exposure]]) || is.character(data[[exposure]]))) {
    stop("exposure must be continuous, binary, or ordinal for PMM.")
  }
  if (anyNA(data[[id_name]])) stop("Subject ids must not be missing.")
  if (!is.numeric(data[[time_name]]) || any(!is.finite(data[[time_name]])) ||
      !is.numeric(data[[response_name]]) || any(!is.finite(data[[response_name]]))) {
    stop("PMM requires finite numeric outcomes and visit times.")
  }

  # Use the sampling recorded in the design, also when raw phase-I data is supplied.
  if (is.null(sampled_col)) {
    sampled_col <- ".mi_sampled"
    design_ids <- as.character(design$data[[id_name]])
    sampled_by_id <- split(design$data$sampled, design_ids)
    if (any(vapply(sampled_by_id, function(x) length(unique(x)) != 1L, logical(1)))) {
      stop("Phase-II sampling must be constant within each subject.")
    }
    sampled_by_id <- vapply(sampled_by_id, function(x) x[1L], numeric(1))
    data[[sampled_col]] <- sampled_by_id[as.character(data[[id_name]])]
  }
  if (!sampled_col %in% names(data)) stop("sampled_col was not found in data.")
  sampled <- data[[sampled_col]]
  if (anyNA(sampled) || any(!sampled %in% c(0, 1))) {
    stop("Phase-II sampling indicators must be 0 or 1 without missing values.")
  }

  data <- data[order(data[[id_name]], data[[time_name]]), , drop = FALSE]
  rows_by_id <- split(seq_len(nrow(data)), as.character(data[[id_name]]))
  times <- sort(unique(data[[time_name]]))
  balanced <- vapply(rows_by_id, function(idx) {
    identical(data[[time_name]][idx], times)
  }, logical(1))
  if (length(times) < 2L || !all(balanced)) {
    stop("PMM requires the same visit times for every subject, without duplicate visits.")
  }
  sampled_by_id <- vapply(rows_by_id, function(idx) {
    x <- data[[sampled_col]][idx]
    if (length(unique(x)) != 1L) stop("Phase-II sampling must be constant within each subject.")
    as.numeric(x[1L])
  }, numeric(1))
  if (!any(sampled_by_id == 1)) stop("No phase-II subjects were identified.")
  if (!any(sampled_by_id == 0)) stop("No phase-I-only subjects were identified for imputation.")
  for (idx in rows_by_id[sampled_by_id == 1]) {
    x <- data[[exposure]][idx]
    if (anyNA(x) || (is.numeric(x) && any(!is.finite(x)))) {
      stop("The expensive covariate has missing or non-finite values among phase-II subjects.")
    }
    if (length(unique(x)) != 1L) stop("The expensive covariate must be constant within each phase-II subject.")
  }

  # One row per subject, with all longitudinal outcomes in the PMM model.
  if (is.null(exposure_model)) {
    exposure_model <- .mi_default_exposure_model(formula, exposure, time_name)
  }
  exposure_rhs <- stats::as.formula(paste("~", .mi_formula_rhs(exposure_model)),
                                    env = environment(exposure_model))
  y_names <- paste0("y", seq_along(times) - 1L)
  exposure_vars <- all.vars(exposure_rhs)
  covariates <- setdiff(exposure_vars, y_names)
  if (any(covariates %in% c(exposure, id_name, time_name, response_name))) {
    stop("The PMM exposure model must use subject-level covariates and wide outcomes y0, y1, ... .")
  }
  if (length(setdiff(covariates, names(data))) > 0L) {
    stop("Missing exposure model column(s): ", paste(setdiff(covariates, names(data)), collapse = ", "))
  }
  if (any(y_names %in% c(exposure, id_name, covariates))) {
    stop("Wide outcome names y0, y1, ... conflict with subject-level column names.")
  }
  for (covariate in covariates) {
    constant <- vapply(rows_by_id, function(idx) {
      length(unique(data[[covariate]][idx])) == 1L
    }, logical(1))
    if (!all(constant)) stop("PMM covariates must be constant within each subject: ", covariate)
  }
  subj_dat <- data[!duplicated(data[[id_name]]),
                   unique(c(id_name, exposure, sampled_col, covariates)), drop = FALSE]
  Ywide <- matrix(data[[response_name]], ncol = length(times), byrow = TRUE)
  colnames(Ywide) <- y_names
  subj_dat <- cbind(subj_dat, Ywide)
  pmm_formula <- stats::update(exposure_rhs,
                              paste("~ . +", paste(y_names, collapse = " + ")))
  pmm_frame <- stats::model.frame(pmm_formula, data = subj_dat, na.action = na.action)
  if (nrow(pmm_frame) != nrow(subj_dat)) stop("PMM predictors must be observed for every subject.")
  H <- stats::model.matrix(pmm_formula, pmm_frame)
  if (any(!is.finite(H))) stop("PMM predictors must be finite.")
  observed <- subj_dat[[sampled_col]] == 1
  X_obs <- subj_dat[[exposure]][observed]
  X_model <- X_obs
  if (is.factor(X_obs)) {
    if (!is.ordered(X_obs) && nlevels(X_obs) != 2L) {
      stop("Use an ordered factor with explicit levels for an ordinal exposure.")
    }
    X_model <- as.numeric(X_obs)
    if (!is.ordered(X_obs)) X_model <- X_model - 1L
  } else if (is.logical(X_obs)) {
    X_model <- as.numeric(X_obs)
  } else if (is.character(X_obs)) {
    binary_levels <- sort(unique(X_obs))
    if (length(binary_levels) != 2L) {
      stop("Character exposures must have two observed categories; use an ordered factor for ordinal exposures.")
    }
    X_model <- match(X_obs, binary_levels) - 1L
  }
  exposure_fit <- .mi_pmm_model(H[observed, , drop = FALSE], X_model,
                               H[!observed, , drop = FALSE], ridge)
  exposure_fit$X_obs <- X_obs
  exposure_fit$formula <- pmm_formula
  exposure_fit$visit_times <- times

  coef_mi <- cov_mi <- completed_data <- NULL
  cov_mi <- vector("list", M)
  if (keep_completed) completed_data <- vector("list", M)
  missing_rows <- data[[sampled_col]] == 0
  recipient_ids <- as.character(subj_dat[[id_name]][!observed])
  fit <- NULL

  for (m in seq_len(M)) {
    if (verbose > 0) message("MI draw ", m, " of ", M, ".")
    exposure_imp <- .mi_pmm_draw(exposure_fit, n_donors)
    imp_map <- setNames(exposure_imp, recipient_ids)
    imp_dat <- data
    imp_dat[[exposure]][missing_rows] <- imp_map[as.character(imp_dat[[id_name]][missing_rows])]
    imp_dat$.mi_weight_all <- 1
    imp_dat$.mi_sampled_all <- 1
    completed_design <- .mi_build_design_like(design, imp_dat, ".mi_weight_all", ".mi_sampled_all")
    fixed_frame <- stats::model.frame(formula, data = imp_dat, na.action = na.action,
                                      drop.unused.levels = TRUE)
    if (nrow(fixed_frame) != nrow(imp_dat)) stop("Completed outcome model covariates must be observed for every visit.")
    fit <- WL(formula, completed_design, init = init, verbose = verbose)
    fit$model.matrix <- stats::model.matrix(formula, fixed_frame)
    if (!fit$Code %in% c(1L, 2L)) stop("Completed-data WL fit did not converge in imputation ", m, ".")
    if (any(!is.finite(fit$coefficients)) || any(!is.finite(fit$covariance)) ||
        any(diag(fit$covariance) <= 0)) stop("Invalid completed-data estimates in imputation ", m, ".")
    if (m == 1L) {
      coef_names <- .mi_coef_names(fit)
      coef_mi <- matrix(NA_real_, M, length(coef_names), dimnames = list(NULL, coef_names))
    }
    coef_mi[m, ] <- coef(fit)
    cov_mi[[m]] <- fit$covariance
    dimnames(cov_mi[[m]]) <- list(coef_names, coef_names)
    if (keep_completed) completed_data[[m]] <- imp_dat
  }
  rubin <- .mi_rubin_multivariate(coef_mi, cov_mi)

  list(
    coefficients   = rubin$beta,
    covariance     = rubin$vcov,
    U_bar          = rubin$U_bar,
    B              = rubin$B,
    coef_mi        = coef_mi,
    cov_mi         = cov_mi,
    exposure_fit   = exposure_fit,
    exposure_model = exposure_model,
    completed_data = completed_data,
    model.matrix   = fit$model.matrix,
    response       = data[[response_name]],
    ids            = data[[id_name]],
    data           = data,
    M              = M,
    n_donors       = n_donors,
    ridge          = ridge
  )
}

#' Fit model using multiple imputation (MI)
#'
#' Impute a continuous, binary, or ordinal subject-level expensive covariate using
#' Type-1 predictive mean matching (PMM). Fit an unweighted ridge linear model among phase-II
#' donors using subject-level covariates and all phase-I longitudinal outcomes.
#' Donor predictions use the fitted coefficients; recipient predictions use
#' posterior draws of the coefficients and residual variance. Each recipient
#' receives the observed exposure of a randomly selected nearest donor.
#' Fit each completed phase-I data set with unit-weight \code{\link{WL}}
#' and combine all model parameters on the unconstrained scale using Rubin's rules.
#' Completed-data sampling probabilities are all 1, so design construction
#' does not introduce additional random sampling between imputation draws.
#'
#' @param formula \code{formula}; a symbolic description of the outcome model,
#'   such as \code{y ~ grp * time + conf}. Use formula interactions to ensure
#'   exposure-dependent terms are recomputed after imputation.
#' @param design an object of class \code{odsdesign} created by \code{\link{ods}}.
#' @param subset an optional logical vector selecting phase-I observations.
#'   For PMM, selected subjects must retain the same complete visit schedule.
#' @param na.action a function indicating what should happen when model
#'   covariates contain NAs. PMM requires observed predictors and outcomes;
#'   missing exposure values are allowed for phase-I-only subjects.
#' @param verbose \code{numeric(1)}; progress printing level.
#' @param init \code{numeric}; initial values passed to completed-data \code{WL}
#'   fits, on the same scale as for \code{WL}.
#' @param ProfileCol reserved for profile likelihood; must be \code{NA} for MI.
#' @param exposure name of the expensive covariate, quoted or unquoted.
#'   Supports numeric continuous, count, binary, or ordinal values, logical
#'   binary values, two-level factors or binary character labels, and ordered
#'   factors with explicitly specified levels for ordinal categories.
#' @param exposure_model a one-sided or two-sided formula specifying
#'   subject-level PMM predictors, e.g. \code{~ conf}. All longitudinal outcomes
#'   are automatically added as \code{y0}, \code{y1}, etc., in increasing order
#'   of visit time. If omitted, covariates come from the outcome formula after
#'   removing exposure and time.
#' @param M number of imputations, an integer of at least 2.
#' @param n_donors number of nearest donors, a positive integer. Defaults to 5;
#'   capped at the number of available donors.
#' @param ridge positive ridge factor for the PMM posterior covariance.
#'   Defaults to \code{1e-4}.
#' @param data optional full phase-I data. Defaults to \code{design$data}.
#' @param sampled_col optional name of a 0/1 phase-II sampling column.
#'   By default PMM uses the sampling recorded in \code{design$data$sampled},
#'   matched by subject id.
#' @param seed optional nonnegative integer random seed.
#' @param keep_completed logical; keep each completed data set if \code{TRUE}.
#' @param ... reserved for additional arguments. Currently no additional
#'   arguments are accepted.
#' @details PMM requires one exposure constant within each donor and
#'   a common complete visit schedule for all subjects. There must be more
#'   donors than coefficients in the PMM working model. The PMM residual sum
#'   of squares must be positive and finite. PMM supports mean,
#'   intercept, slope, and bivariate ODS designs with \code{method},
#'   \code{p_sample}, and cutpoints stored in the design. Imputed values stay
#'   on the observed donor support and retain the original exposure type,
#'   factor levels, and ordering. Numeric exposures retain their numeric scores
#'   in the working model; logical and two-level factor exposures use 0/1,
#'   binary character labels use 0/1 in sorted label order, and ordered factors
#'   use scores 1, 2, etc. in the declared level order. These numeric scores are
#'   only used for matching; recipients receive the donor's original value.
#'   Neither sampling weights nor a Poisson regression are used in the PMM
#'   working model. Completed-data analyses use model-based WL covariance.
#' @return An object of class \code{"MI"} containing pooled
#'   \code{coefficients}, \code{covariance}, within-imputation covariance
#'   \code{U_bar}, between-imputation covariance \code{B}, individual
#'   \code{coef_mi} and \code{cov_mi}, the PMM \code{exposure_fit}, and optional
#'   \code{completed_data}. The methods \code{coef}, \code{fixef}, \code{ranef},
#'   \code{vcov}, \code{print}, and \code{summary} follow the fitted-model interface.
#' @export
#' @examples
#' \dontrun{
#' fit <- MI(y ~ grp * time + conf, design,
#'           exposure = grp, exposure_model = ~ conf,
#'           M = 20, n_donors = 5, seed = 1)
#' coef(fit)
#' summary(fit)
#' }
MI <- function(
    formula,
    design,
    subset = NULL,
    na.action = getOption('na.action'),
    verbose = 0L,
    init = NULL,
    ProfileCol = NA,
    exposure,
    exposure_model = NULL,
    M = 20L,
    n_donors = 5L,
    ridge = 1e-4,
    data = NULL,
    sampled_col = NULL,
    seed = NULL,
    keep_completed = FALSE,
    ...)
{
  # Validate arguments
  coll <- makeAssertCollection()
  assert_formula(formula, add=coll)
  assert_class(design, "odsdesign", add=coll)
  assert_numeric(M, len=1L, lower=2, finite=TRUE, any.missing=FALSE, add=coll)
  assert_numeric(n_donors, len=1L, lower=1, finite=TRUE, any.missing=FALSE, add=coll)
  assert_numeric(ridge, len=1L, lower=.Machine$double.eps, finite=TRUE, any.missing=FALSE, add=coll)
  assert_logical(keep_completed, len=1L, any.missing=FALSE, add=coll)
  assert_formula(exposure_model, null.ok=TRUE, add=coll)
  assert_numeric(seed, len=1L, lower=0, upper=.Machine$integer.max,
                 finite=TRUE, any.missing=FALSE, null.ok=TRUE, add=coll)
  reportAssertions(coll)
  if (length(formula) != 3L) stop("formula must be two-sided.")
  if (M != floor(M)) stop("M must be an integer of at least 2.")
  if (n_donors != floor(n_donors)) stop("n_donors must be a positive integer.")
  if (!is.null(seed) && seed != floor(seed)) stop("seed must be an integer.")
  if (length(ProfileCol) != 1L || !is.na(ProfileCol) || any(!is.na(design$ProfileCol))) {
    stop("ProfileCol must be NA for MI.")
  }
  if (length(verbose) != 1L || is.na(verbose) ||
      !(is.numeric(verbose) || is.logical(verbose)) || !is.finite(verbose) || verbose < 0) {
    stop("verbose must be a nonnegative printing level.")
  }
  if (missing(exposure)) stop("Please provide the expensive covariate name in exposure.")
  exposure_expr <- substitute(exposure)
  exposure_name <- if (is.name(exposure_expr)) as.character(exposure_expr) else exposure
  assert_character(exposure_name, len=1L, any.missing=FALSE)
  assert_character(sampled_col, len=1L, any.missing=FALSE, null.ok=TRUE)
  if (is.null(data)) data <- design$data
  if (!is.data.frame(data) || nrow(data) == 0L) stop("data must be a nonempty phase-I data frame.")
  data <- as.data.frame(data)

  cl <- match.call()
  subset <- eval(substitute(subset), data, parent.frame())
  if (!is.null(subset)) {
    assert_logical(subset, len=nrow(data), any.missing=FALSE)
    data <- data[subset, , drop = FALSE]
    if (nrow(data) == 0L) stop("subset must retain phase-I observations.")
  }
  if (length(list(...)) > 0L) stop("Unused additional arguments supplied to MI.")
  if (!is.null(seed)) set.seed(seed)
  fit <- MI_internal(
    formula = formula, design = design, exposure = exposure_name,
    exposure_model = exposure_model, M = as.integer(M), n_donors = as.integer(n_donors),
    ridge = ridge, data = data, sampled_col = sampled_col, na.action = na.action,
    verbose = verbose, init = init, keep_completed = keep_completed
  )
  fit$formula <- formula
  fit$design <- design
  fit$exposure <- exposure_name
  fit$call <- cl
  class(fit) <- "MI"
  fit
}

#' @export
fixef.MI <- function(object, ...)
{
  est <- object$coefficients
  est[seq_len(length(est) - triangle(object$design$n_rand) - 1L)]
}

#' @export
ranef.MI <- function(object, transform = FALSE, ...)
{
  est <- object$coefficients
  le <- length(est)
  ran <- est[(le - triangle(object$design$n_rand)):le]

  if (transform) ranef_transform(ran, object$design$n_rand) else ran
}

#' @exportS3Method
coef.MI <- function(object, complete = TRUE, transform = FALSE, ...)
{
  c(fixef(object, ...), ranef(object, transform = transform, ...))
}

#' @exportS3Method
vcov.MI <- function(object, complete = TRUE, ...)
{
  nm <- names(coef(object))
  vc <- object$covariance
  rownames(vc) <- nm
  colnames(vc) <- nm
  vc
}

#' @exportS3Method
print.MI <- function(x, digits = max(3L, getOption("digits")), transform = FALSE, ...)
{
  object <- x
  cat("\nCall:\n",
      paste(deparse(object$call), collapse = "\n"),
      "\n\n",
      "Expensive covariate: ", object$exposure,
      "\n",
      "MI draws: ", object$M,
      "\n",
      sep = "")
  cat("\nFixed Effects:\n")
  print(round(fixef(object), digits = digits), ...)
  cat("\nRandom Effects:\n")
  print(round(ranef(object, transform = transform), digits = digits), ...)
  invisible(object)
}

#' @exportS3Method
summary.MI <- function(object, digits = max(3L, getOption("digits")),
                       transform = FALSE, ...)
{
  raw <- coef(object, ...)
  beta <- coef(object, transform = transform, ...)
  se <- sqrt(diag(vcov(object)))
  object$transform <- transform
  object$digits <- digits
  z <- raw / se
  le <- length(beta)
  object$n_random <- triangle(object$design$n_rand) + 1L
  object$n_fixed <- le - object$n_random

  object$coefficients <- cbind(
    Estimate = beta,
    `Std. Error` = se,
    L95 = raw + se * stats::qnorm(0.025),
    U95 = raw + se * stats::qnorm(0.975),
    `z value` = z,
    `Pr(>|z|)` = 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  )

  if (transform) {
    n_rand <- object$design$n_rand
    ran <- (le - triangle(n_rand)):nrow(object$coefficients)
    object$coefficients[ran, "L95"] <- ranef_transform(object$coefficients[ran, "L95"], n_rand)
    object$coefficients[ran, "U95"] <- ranef_transform(object$coefficients[ran, "U95"], n_rand)

    rhos <- ran[(n_rand + 1L):(length(ran) - 1L)]
    object$coefficients[rhos, 2L] <- abs(1 / cosh(raw[rhos] / 2)^2 / 2) *
      object$coefficients[rhos, 2L]

    ran <- ran[!ran %in% rhos]
    object$coefficients[ran, 2L] <- object$coefficients[ran, 1L] *
      object$coefficients[ran, 2L]
  }

  class(object) <- c("summary.MI", "MI")
  object
}

#' @exportS3Method
print.summary.MI <- function(x, digits = NULL,
                             signif.stars = FALSE, ...)
{
  object <- x
  if (is.null(digits)) digits <- object$digits

  cat("\nCall:\n",
      paste(deparse(object$call), collapse = "\n"),
      "\n\n",
      "Expensive covariate: ", object$exposure,
      "\n",
      "MI draws: ", object$M,
      "\n",
      sep = "")
  cat("\nFixed Effects:\n")
  stats::printCoefmat(
    object$coefficients[seq_len(object$n_fixed), , drop = FALSE],
    digits = digits - 1L,
    dig.tst = digits,
    signif.stars = signif.stars,
    na.print = "NA",
    ...
  )

  cat("\nRandom Effects:\n")
  random <- object$coefficients[(object$n_fixed + 1L):nrow(object$coefficients), , drop = FALSE]
  nt <- max(nchar(object$design$id) + 1L, 7L)
  pad <- paste0(rep(" ", nt), collapse = "")
  rownames(random) <-
    if (object$transform) {
      c(
        paste(object$design$id, "(Intercept)"),
        paste0(pad, object$design$time),
        paste0(pad, "rho"),
        "Residual"
      )
    } else {
      c(
        paste(object$design$id, "log(Intercept)"),
        paste0(pad, "log(", object$design$time, ")"),
        paste0(pad, "2*atanh(rho)"),
        "log(Residual)"
      )
    }

  print(round(random[, 1:4, drop = FALSE], digits))
  invisible(object)
}
