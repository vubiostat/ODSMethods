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

.mi_weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]

  if (length(x) < 2) stop("Need at least two finite observations for weighted quantiles.")
  if (any(probs < 0 | probs > 1)) stop("probs must be in [0, 1].")

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)

  stats::approx(
    x = c(0, cw),
    y = c(x[1], x),
    xout = probs,
    method = "constant",
    f = 0,
    rule = 2,
    ties = "ordered"
  )$y
}

.mi_make_intervals <- function(x_phase2, w_phase2, L = 20L) {
  ok <- is.finite(x_phase2) & is.finite(w_phase2) & w_phase2 > 0
  x_phase2 <- x_phase2[ok]
  w_phase2 <- w_phase2[ok]

  if (L < 3) stop("L must be at least 3.")
  if (length(x_phase2) < L) {
    stop("Not enough phase-II subjects to create ", L, " weighted-quantile intervals.")
  }

  probs <- seq(0, 1, length.out = L + 1L)
  breaks0 <- as.numeric(.mi_weighted_quantile(x_phase2, w_phase2, probs))

  if (any(!is.finite(breaks0))) stop("Weighted quantile cutpoints contain non-finite values.")
  if (any(diff(breaks0) <= 0)) {
    stop(
      "Weighted quantile cutpoints are not strictly increasing. ",
      "Try a smaller L or check whether the observed phase-II exposure has many ties."
    )
  }

  breaks <- breaks0
  breaks[1L] <- -Inf
  breaks[length(breaks)] <- Inf

  interval <- cut(
    x_phase2,
    breaks = breaks,
    labels = seq_len(L),
    include.lowest = TRUE,
    right = FALSE
  )
  interval <- as.integer(as.character(interval))

  tab <- table(factor(interval, levels = seq_len(L)))
  if (any(tab == 0)) {
    stop(
      "At least one weighted-quantile interval is empty in phase-II data. ",
      "Try a smaller L."
    )
  }

  midpoints <- rep(NA_real_, L)
  for (ll in seq_len(L)) {
    idx <- which(interval == ll)
    midpoints[ll] <- sum(w_phase2[idx] * x_phase2[idx]) / sum(w_phase2[idx])
  }

  list(
    L = L,
    breaks = breaks,
    cutpoints_finite = breaks0,
    midpoints = midpoints,
    levels = seq_len(L),
    phase2_counts = as.integer(tab)
  )
}

.mi_assign_interval <- function(x, interval_info) {
  out <- cut(
    x,
    breaks = interval_info$breaks,
    labels = interval_info$levels,
    include.lowest = TRUE,
    right = FALSE
  )
  as.integer(as.character(out))
}

.mi_get_covariance <- function(fit) {
  if (is.null(fit$covariance)) {
    stop("The fitted WL object must contain a covariance component.")
  }
  fit$covariance
}

.mi_get_orm_vcov_all <- function(fit) {
  tryCatch(
    stats::vcov(fit, intercepts = "all"),
    error = function(e) stats::vcov(fit)
  )
}

.mi_draw_mvn <- function(mu, Sigma) {
  Sigma <- as.matrix(Sigma)
  if (length(mu) != nrow(Sigma) || length(mu) != ncol(Sigma)) {
    stop("Coefficient vector and covariance matrix dimensions do not match.")
  }

  if (is.null(names(mu))) names(mu) <- paste0("par", seq_along(mu))
  if (is.null(rownames(Sigma))) rownames(Sigma) <- names(mu)
  if (is.null(colnames(Sigma))) colnames(Sigma) <- names(mu)

  Sigma <- Sigma[names(mu), names(mu), drop = FALSE]
  Sigma <- as.matrix(Matrix::nearPD(Sigma)$mat)
  out <- as.numeric(MASS::mvrnorm(1, mu = mu, Sigma = Sigma))
  names(out) <- names(mu)
  out
}

.mi_get_orm_probs <- function(fit, newdata, interval_levels) {
  p <- stats::predict(fit, newdata = newdata, type = "fitted.ind")
  p <- as.numeric(p)

  if (length(p) != length(interval_levels)) {
    stop(
      "The number of predicted exposure probabilities from orm is ", length(p),
      ", but the number of candidate intervals is ", length(interval_levels), "."
    )
  }

  p <- pmax(p, .Machine$double.xmin)
  p <- p / sum(p)
  names(p) <- as.character(interval_levels)
  p
}

.mi_rubin_multivariate <- function(coef_mi, cov_mi) {
  M <- nrow(coef_mi)
  beta_bar <- colMeans(coef_mi)
  U_bar <- Reduce("+", cov_mi) / M

  B <- matrix(0, ncol(coef_mi), ncol(coef_mi))
  for (mm in seq_len(M)) {
    d <- coef_mi[mm, ] - beta_bar
    B <- B + tcrossprod(d)
  }
  B <- B / (M - 1)

  total_var <- U_bar + (1 + 1 / M) * B
  dimnames(total_var) <- list(colnames(coef_mi), colnames(coef_mi))
  dimnames(U_bar) <- list(colnames(coef_mi), colnames(coef_mi))
  dimnames(B) <- list(colnames(coef_mi), colnames(coef_mi))

  list(beta = beta_bar, vcov = total_var, U_bar = U_bar, B = B)
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
    p_sample = design$p_sample,
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

#' Multiple imputation analysis for ODS data
#'
#' Fit a multiple-imputation analysis for a continuous expensive covariate in
#' an outcome dependent sampling design. The observed phase-II expensive
#' covariate is discretized into weighted quantile intervals, an ordinal
#' exposure model is fit with \code{rms::orm}, unsampled phase-I subjects are
#' imputed, and the completed-data analyses are combined with Rubin's rule.
#'
#' @param formula A model formula for the outcome model, such as
#'   \code{y ~ grp * time + conf}.
#' @param design An \code{odsdesign} object created by \code{\link{ods}}.
#' @param exposure The name of the continuous expensive covariate to impute.
#'   This can be supplied quoted or unquoted.
#' @param exposure_model A one-sided formula for the subject-level exposure
#'   model, such as \code{~ conf}. If omitted, subject-level covariates are
#'   taken from the outcome model after removing the exposure and time
#'   variables.
#' @param M Number of final MI draws.
#' @param n_imp Number of updating imputations within each final draw.
#' @param L Number of weighted quantile intervals used to approximate the
#'   continuous expensive covariate.
#' @param data Optional full phase-I data. Defaults to \code{design$data}.
#' @param sampled_col Name of the column indicating phase-II sampling. Defaults
#'   to \code{"sampled"}.
#' @param seed Optional random seed.
#' @param verbose Logical; print progress messages.
#' @param init Optional initial values passed to the initial \code{\link{WL}}
#'   fit.
#' @param keep_completed Logical; if \code{TRUE}, keep the completed data set
#'   from each final MI draw.
#' @param ... Reserved for future options.
#'
#' @return An object of class \code{"MI"} with combined coefficients,
#'   covariance matrix, the initial WL fit, the ordinal exposure fit, and
#'   interval information.
#' @export
#'
#' @examples
#' \dontrun{
#' fit_mi <- MI(
#'   y ~ grp * time + conf,
#'   design,
#'   exposure = grp,
#'   exposure_model = ~ conf,
#'   M = 20,
#'   n_imp = 5,
#'   L = 20
#' )
#' summary(fit_mi)
#' }
MI <- function(
    formula,
    design,
    exposure,
    exposure_model = NULL,
    M = 20L,
    n_imp = 5L,
    L = 20L,
    data = NULL,
    sampled_col = NULL,
    seed = NULL,
    verbose = FALSE,
    init = NULL,
    keep_completed = FALSE,
    ...)
{
  if (!inherits(design, "odsdesign")) stop("design must be an odsdesign object.")
  if (missing(exposure)) stop("Please provide the expensive covariate name in exposure.")

  exposure_expr <- substitute(exposure)
  exposure_name <- if (is.character(exposure_expr)) exposure_expr[1L] else as.character(exposure_expr)
  weights_name <- design$weights

  M <- as.integer(M)
  n_imp <- as.integer(n_imp)
  L <- as.integer(L)

  if (length(M) != 1L || M < 2L) stop("M must be at least 2.")
  if (length(n_imp) != 1L || n_imp < 1L) stop("n_imp must be at least 1.")
  if (length(L) != 1L || L < 3L) stop("L must be at least 3.")

  if (!is.null(seed)) set.seed(seed)

  if (is.null(data)) data <- design$data
  data <- as.data.frame(data)

  id_name <- design$id
  time_name <- design$time
  response_name <- design$response

  needed <- c(id_name, time_name, response_name, exposure_name)
  missing_needed <- setdiff(needed, names(data))
  if (length(missing_needed) > 0L) {
    stop("Missing required column(s): ", paste(missing_needed, collapse = ", "))
  }

  if (is.null(sampled_col)) sampled_col <- "sampled"
  if (!sampled_col %in% names(data)) {
    stop("sampled_col '", sampled_col, "' was not found in data.")
  }

  if (is.null(weights_name)) {
    weights_name <- ".mi_unit_weight"
    data[[weights_name]] <- 1
  } else if (!weights_name %in% names(data)) {
    stop("The weight column stored in design$weights was not found in data: ", weights_name)
  }

  data[[weights_name]] <- as.numeric(data[[weights_name]])
  if (any(!is.finite(data[[weights_name]]) | data[[weights_name]] <= 0, na.rm = TRUE)) {
    stop("weights must be finite and positive.")
  }

  sampled <- data[[sampled_col]] == 1
  sampled[is.na(sampled)] <- FALSE
  id_chr <- as.character(data[[id_name]])
  sampled_by_id <- tapply(sampled, id_chr, any)
  phase2_ids <- names(sampled_by_id)[sampled_by_id]

  if (length(phase2_ids) == 0L) stop("No phase-II subjects were identified.")
  if (length(phase2_ids) == length(sampled_by_id)) {
    stop("No phase-I-only subjects were identified for imputation.")
  }

  phase2_dat <- data[id_chr %in% phase2_ids, , drop = FALSE]
  phase1_only <- data[!id_chr %in% phase2_ids, , drop = FALSE]

  if (any(is.na(phase2_dat[[exposure_name]]))) {
    stop("The expensive covariate has missing values among phase-II subjects.")
  }

  first_phase2 <- !duplicated(as.character(phase2_dat[[id_name]]))
  phase2_subj_for_bins <- phase2_dat[first_phase2, , drop = FALSE]

  interval_info <- .mi_make_intervals(
    x_phase2 = phase2_subj_for_bins[[exposure_name]],
    w_phase2 = phase2_subj_for_bins[[weights_name]],
    L = L
  )
  interval_levels <- interval_info$levels

  phase2_dat[[".mi_xe_int"]] <- .mi_assign_interval(phase2_dat[[exposure_name]], interval_info)
  phase2_subj <- phase2_dat[!duplicated(as.character(phase2_dat[[id_name]])), , drop = FALSE]

  tab_phase2 <- table(factor(phase2_subj[[".mi_xe_int"]], levels = interval_levels))
  if (any(tab_phase2 == 0)) stop("Some exposure intervals are empty in phase-II data.")

  phase2_subj[[".mi_xe_ord"]] <- ordered(phase2_subj[[".mi_xe_int"]], levels = interval_levels)

  if (is.null(exposure_model)) {
    exposure_model <- .mi_default_exposure_model(formula, exposure_name, time_name)
  }
  exposure_rhs <- .mi_formula_rhs(exposure_model)
  exposure_formula <- stats::as.formula(paste(".mi_xe_ord ~", exposure_rhs))

  exposure_vars <- all.vars(exposure_formula)
  missing_exposure_vars <- setdiff(exposure_vars, names(phase2_subj))
  if (length(missing_exposure_vars) > 0L) {
    stop("Missing exposure model column(s): ", paste(missing_exposure_vars, collapse = ", "))
  }

  phase2_fit_dat <- phase2_dat
  phase2_fit_dat[[".mi_sampled_fit"]] <- 1

  if (verbose) message("Fitting initial WL outcome model.")
  phase2_design <- .mi_build_design_like(
    design = design,
    data = phase2_fit_dat,
    weights_name = weights_name,
    sampled_name = ".mi_sampled_fit"
  )
  fit_initial <- WL(formula, phase2_design, init = init)

  if (verbose) message("Fitting ordinal exposure model.")
  fit_exposure <- rms::orm(
    formula = exposure_formula,
    data = phase2_subj,
    weights = phase2_subj[[weights_name]],
    family = "probit",
    x = TRUE,
    y = TRUE
  )

  omega_hat <- stats::coef(fit_exposure)
  V_omega <- .mi_get_orm_vcov_all(fit_exposure)
  if (length(omega_hat) != ncol(V_omega)) {
    stop("Exposure-model coefficient and covariance dimensions do not match.")
  }

  coef_names <- .mi_coef_names(fit_initial)
  if (length(coef_names) != length(fit_initial$coefficients)) {
    stop("Could not construct coefficient names for the MI fit.")
  }
  coef_mi <- matrix(NA_real_, nrow = M, ncol = length(coef_names))
  colnames(coef_mi) <- coef_names
  cov_mi <- vector("list", M)
  completed_data <- if (keep_completed) vector("list", M) else NULL

  outcome_terms <- stats::terms(formula, data = phase2_fit_dat)
  random_vars <- all.vars(design$formula[[3L]][[2L]])
  if (length(random_vars) != 1L) {
    stop("MI currently expects one random time variable in the ods design.")
  }
  if (!identical(random_vars, time_name)) {
    stop("MI currently expects the ods random time variable to match design$time.")
  }

  for (mm in seq_len(M)) {
    if (verbose) message("MI draw ", mm, " of ", M, ".")
    fit_current <- fit_initial
    imp_dat <- NULL

    for (kk in seq_len(n_imp)) {
      params_draw <- .mi_draw_mvn(
        mu = fit_current$coefficients,
        Sigma = fit_current$covariance
      )

      X1 <- fit_current$model.matrix
      p_fix <- ncol(X1)
      beta_all <- params_draw[seq_len(p_fix)]
      names(beta_all) <- colnames(X1)

      if (design$n_rand != 2L) stop("MI currently expects random intercept and random slope.")
      k0 <- p_fix
      sigma_0 <- exp(params_draw[k0 + 1L])
      sigma_1 <- exp(params_draw[k0 + 2L])
      rho <- (exp(params_draw[k0 + 3L]) - 1) / (exp(params_draw[k0 + 3L]) + 1)
      D_draw <- matrix(
        c(
          sigma_0^2, rho * sigma_0 * sigma_1,
          rho * sigma_0 * sigma_1, sigma_1^2
        ),
        nrow = 2L
      )
      sigmae_2 <- exp(params_draw[k0 + 4L])^2

      omega_draw <- .mi_draw_mvn(omega_hat, V_omega)
      fit_exposure_draw <- fit_exposure
      fit_exposure_draw$coefficients <- omega_draw

      phase1_ids <- as.character(phase1_only[[id_name]])
      uid <- unique(phase1_ids)
      exposure_by_id <- rep(NA_real_, length(uid))
      names(exposure_by_id) <- uid

      for (ii in seq_along(uid)) {
        id_i <- uid[ii]
        idx <- which(phase1_ids == id_i)
        dat_i <- phase1_only[idx, , drop = FALSE]

        Zi <- cbind(1, dat_i[[time_name]])
        Vi <- Zi %*% D_draw %*% t(Zi) + sigmae_2 * diag(nrow(Zi))
        y_i <- matrix(dat_i[[response_name]], ncol = 1L)

        new_subj <- dat_i[1L, , drop = FALSE]
        px_i <- .mi_get_orm_probs(
          fit = fit_exposure_draw,
          newdata = new_subj,
          interval_levels = interval_levels
        )

        log_w <- rep(NA_real_, length(interval_levels))
        names(log_w) <- as.character(interval_levels)

        for (ll in seq_along(interval_levels)) {
          interval_val <- interval_levels[ll]
          x_val <- interval_info$midpoints[interval_val]
          dat_l <- dat_i
          dat_l[[exposure_name]] <- x_val

          X_l <- stats::model.matrix(outcome_terms, data = dat_l)
          missing_cols <- setdiff(names(beta_all), colnames(X_l))
          if (length(missing_cols) > 0L) {
            stop("Imputed model matrix is missing column(s): ", paste(missing_cols, collapse = ", "))
          }
          X_l <- X_l[, names(beta_all), drop = FALSE]

          mu_l <- matrix(X_l %*% beta_all, ncol = 1L)
          res_l <- y_i - mu_l

          log_det <- as.numeric(determinant(Vi, logarithm = TRUE)$modulus)
          quad <- as.numeric(t(res_l) %*% solve(Vi, res_l))
          log_f_y <- -0.5 * (length(y_i) * log(2 * pi) + log_det + quad)

          log_w[ll] <- log_f_y + log(px_i[as.character(interval_val)])
        }

        prob_l <- exp(log_w - max(log_w))
        prob_l <- prob_l / sum(prob_l)
        if (any(!is.finite(prob_l)) || sum(prob_l) <= 0) {
          stop("Could not compute finite imputation probabilities for subject ", id_i, ".")
        }

        interval_draw <- as.integer(sample(names(prob_l), size = 1L, prob = as.numeric(prob_l)))
        exposure_by_id[id_i] <- interval_info$midpoints[interval_draw]
      }

      phase1_only_imp <- phase1_only
      phase1_only_imp[[exposure_name]] <- as.numeric(exposure_by_id[as.character(phase1_only_imp[[id_name]])])

      phase2_completed <- phase2_dat
      phase2_completed[[".mi_xe_int"]] <- NULL
      phase1_only_imp[[".mi_sampled_all"]] <- 1
      phase2_completed[[".mi_sampled_all"]] <- 1
      phase1_only_imp[[".mi_weight_all"]] <- 1
      phase2_completed[[".mi_weight_all"]] <- 1

      imp_dat <- rbind(phase1_only_imp, phase2_completed)
      imp_dat <- imp_dat[order(imp_dat[[id_name]], imp_dat[[time_name]]), , drop = FALSE]

      completed_design <- .mi_build_design_like(
        design = design,
        data = imp_dat,
        weights_name = ".mi_weight_all",
        sampled_name = ".mi_sampled_all"
      )
      fit_current <- WL(formula, completed_design)
    }

    coef_mi[mm, ] <- fit_current$coefficients
    cov_mi[[mm]] <- .mi_get_covariance(fit_current)
    dimnames(cov_mi[[mm]]) <- list(coef_names, coef_names)
    if (keep_completed) completed_data[[mm]] <- imp_dat
  }

  rubin <- .mi_rubin_multivariate(coef_mi, cov_mi)
  names(rubin$beta) <- coef_names

  out <- list(
    call = match.call(),
    formula = formula,
    design = design,
    exposure = exposure_name,
    exposure_model = exposure_model,
    coefficients = rubin$beta,
    covariance = rubin$vcov,
    U_bar = rubin$U_bar,
    B = rubin$B,
    coef_mi = coef_mi,
    cov_mi = cov_mi,
    interval_info = interval_info,
    initial_fit = fit_initial,
    exposure_fit = fit_exposure,
    completed_data = completed_data,
    M = M,
    n_imp = n_imp,
    L = L
  )
  class(out) <- "MI"
  out
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
    object$coefficients[rhos, 2L] <- abs(1 / cosh(coef(object)[rhos] / 2)^2 / 2) *
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
                             signif.stars = getOption("show.signif.stars"), ...)
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
