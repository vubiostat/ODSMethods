# PMM algorithm reference from optimal_vs_bds_pmm_type1_yanb2.R. The original data comes from ordmi_slp_pois_pmm_type1_yanb1.R;

validated_mi_data <- function(N1 = 160L, seed = 20261001, continuous = FALSE,
                              exposure_type = NULL)
{
  set.seed(seed)
  n <- 6L
  conf <- rnorm(N1)
  grp <- if (identical(exposure_type, "binary")) {
    rbinom(N1, 1, plogis(0.35 * conf))
  } else if (identical(exposure_type, "ordinal")) {
    as.integer(cut(rnorm(N1, mean = 0.35 * conf), c(-Inf, -0.5, 0.5, Inf)))
  } else if (continuous) rnorm(N1, mean = 1 + 0.35 * conf, sd = 1) else
    rpois(N1, lambda = exp(1 + 0.35 * conf))
  b_i <- MASS::mvrnorm(N1, mu = c(0, 0), Sigma = diag(c(9^2, 1.25^2)))
  dat <- data.frame(id = rep(seq_len(N1), each = n),
                    time = rep(0:(n - 1L), times = N1),
                    grp = rep(grp, each = n), conf = rep(conf, each = n))
  dat$y <- 75 - dat$grp - 0.5 * dat$time - 2 * dat$conf -
    0.25 * dat$grp * dat$time + rep(b_i[, 1], each = n) +
    rep(b_i[, 2], each = n) * dat$time + rnorm(nrow(dat), sd = 3.5)
  dat$WeightsSlp <- 1
  dat
}

validated_mi_design <- function(method = "slope", data = validated_mi_data())
{
  force(data)
  set.seed(20261001)
  cutpoints <- if (method == "bivariate") c(60, 85, -3.2561525, 0.7995312) else
    if (method == "slope") c(-3.2561525, 0.7995312) else c(60, 85)
  design <- ods(y ~ time|id, method,
                p_sample = if (method == "bivariate") c(0.125, 1) else c(1, 0.125, 1),
                data = data, cutpoints = cutpoints, weights = "WeightsSlp")
  if (method == "bivariate") {
    inside <- design$z_i["intercept", ] > cutpoints[1] &
      design$z_i["intercept", ] < cutpoints[2] &
      design$z_i["slope", ] > cutpoints[3] & design$z_i["slope", ] < cutpoints[4]
    prob <- ifelse(inside, 0.125, 1)
  } else {
    prob <- design$p_sample[as.integer(cut(design$z_i[method, ],
                                           c(-Inf, cutpoints, Inf)))]
  }
  names(prob) <- colnames(design$z_i)
  design$data$WeightsSlp <- 1 / prob[as.character(design$data$id)]
  design
}

validated_mi <- function(design, M = 3L, n_donors = 5L, ridge = 1e-4, seed = 37L)
{
  force(design)
  set.seed(seed)
  phase1_dat <- design$data[order(design$data$id, design$data$time), ]
  phase1_dat$SampledDesign <- phase1_dat$sampled
  n_visit <- length(unique(phase1_dat$time))
  subj_dat <- phase1_dat[!duplicated(phase1_dat$id),
                         c("id", "grp", "conf", "SampledDesign")]
  Ywide <- matrix(phase1_dat$y, ncol = n_visit, byrow = TRUE)
  y_names <- paste0("y", 0:(n_visit - 1))
  colnames(Ywide) <- y_names
  subj_dat <- cbind(subj_dat, Ywide)
  donor <- subj_dat[subj_dat$SampledDesign == 1, , drop = FALSE]
  recipient <- subj_dat[subj_dat$SampledDesign == 0, , drop = FALSE]
  if (nrow(donor) <= n_visit + 2) stop("Too few Phase-II donors for PMM.")

  pmm_formula <- reformulate(c("conf", y_names))
  H_obs <- model.matrix(pmm_formula, data = donor)
  X_obs <- donor$grp
  S <- crossprod(H_obs)
  V <- solve(S + diag(diag(S) * ridge))
  alpha_hat <- as.numeric(V %*% crossprod(H_obs, X_obs))
  names(alpha_hat) <- colnames(H_obs)
  pred_donor <- as.numeric(H_obs %*% alpha_hat)
  RSS <- sum((X_obs - pred_donor)^2)
  df_resid <- nrow(H_obs) - ncol(H_obs)
  if (!is.finite(RSS) || RSS <= 0 || df_resid <= 0) {
    stop("Invalid PMM RSS or residual degrees of freedom.")
  }
  H_mis <- model.matrix(pmm_formula, data = recipient)
  V_half <- chol(V)
  coef_mi <- matrix(NA_real_, M, 9L)
  cov_mi <- completed_data <- vector("list", M)
  codes <- integer(M)

  for (m in seq_len(M)) {
    g_draw <- rchisq(1, df = df_resid)
    sigma2_draw <- RSS / g_draw
    z1 <- rnorm(length(alpha_hat))
    alpha_draw <- alpha_hat + sqrt(sigma2_draw) * as.numeric(z1 %*% V_half)
    pred_recipient <- as.numeric(H_mis %*% alpha_draw)
    grp_imp <- numeric(nrow(recipient))
    d_use <- min(n_donors, nrow(donor))
    for (i in seq_len(nrow(recipient))) {
      distance_i <- abs(pred_donor - pred_recipient[i])
      nearest_i <- order(distance_i)[seq_len(d_use)]
      selected_donor <- sample(nearest_i, size = 1)
      grp_imp[i] <- donor$grp[selected_donor]
    }
    phase1_dat_imp <- phase1_dat
    imp_map <- setNames(grp_imp, recipient$id)
    miss <- phase1_dat_imp$SampledDesign == 0
    phase1_dat_imp$grp[miss] <- imp_map[as.character(phase1_dat_imp$id[miss])]
    phase1_dat_imp$grp_time <- phase1_dat_imp$grp * phase1_dat_imp$time
    phase1_dat_imp$wt1 <- 1
    phase1_dat_imp$sampled_for_fit <- 1
    design_completed <- ods(y ~ time|id, "slope", p_sample = c(1, 1, 1),
                            data = phase1_dat_imp, quantiles = c(0.1, 0.9),
                            weights = "wt1", sampled_name = "sampled_for_fit")
    fit_m <- WL(y ~ grp * time + conf, design_completed, data = phase1_dat_imp)
    coef_mi[m, ] <- as.numeric(coef(fit_m))
    cov_mi[[m]] <- if (!is.null(fit_m$covariance)) fit_m$covariance else vcov(fit_m)
    codes[m] <- fit_m$Code
    completed_data[[m]] <- phase1_dat_imp
  }
  beta_bar <- colMeans(coef_mi)
  U_bar <- Reduce("+", cov_mi) / M
  B <- cov(coef_mi)
  T_mi <- U_bar + (1 + 1 / M) * B
  list(coefficients = beta_bar, covariance = T_mi,
       coef_mi = coef_mi, cov_mi = cov_mi, U_bar = U_bar, B = B,
       completed_data = completed_data, codes = codes)
}
