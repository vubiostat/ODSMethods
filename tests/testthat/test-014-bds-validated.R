context("BDS Validated Likelihood Methods")

bds_test_reference_ids <- function(data, n_each = 10)
{
  genotype_by_id <- tapply(
    data$Genotype,
    data$Patient,
    function(x) unique(stats::na.omit(x))[1]
  )

  false_ids <- names(genotype_by_id)[!is.na(genotype_by_id) & !genotype_by_id]
  true_ids  <- names(genotype_by_id)[!is.na(genotype_by_id) & genotype_by_id]

  c(head(false_ids, n_each), head(true_ids, n_each))
}

bds_test_fit_data <- function(data, n_each = 15)
{
  ids <- bds_test_reference_ids(data, n_each = n_each)
  dat <- data[as.character(data$Patient) %in% ids, ]
  dat$wt_ref <- 1
  dat
}

bds_test_acml_input <- function(design, ids)
{
  dat <- design$data[as.character(design$data[[design$id]]) %in% ids, ]
  mf <- stats::model.frame(Response ~ Month * Genotype + Patient, data = dat)

  lookup <- design$data[
    !duplicated(as.character(design$data[[design$id]])),
    c(design$id, "cutpoints_i", "acml_samp_prob_i"),
    drop = FALSE
  ]
  rownames(lookup) <- as.character(lookup[[design$id]])

  list(
    y = mf[[design$response]],
    x = stats::model.matrix(Response ~ Month * Genotype, data = mf),
    z = cbind(`(Intercept)` = 1, Month = mf[[design$time]]),
    id = mf[[design$id]],
    cutpoints = do.call(rbind, lookup[as.character(mf[[design$id]]), "cutpoints_i"]),
    SampProb = do.call(rbind, lookup[as.character(mf[[design$id]]), "acml_samp_prob_i"])
  )
}

bds_test_params <- function(x)
{
  beta <- c(
    `(Intercept)` = 10,
    Month = 0.5,
    GenotypeTRUE = -0.5,
    `Month:GenotypeTRUE` = 0.25
  )

  c(
    beta[colnames(x)],
    log(2),
    log(0.5),
    log((1 + 0.1) / (1 - 0.1)),
    log(1)
  )
}

test_that("bds constructs BLUP sampling variables from the phase-I mixed model",
{
  data(gbti)
  set.seed(202405)

  design <- bds(
    Response ~ Month | Patient,
    "slope",
    p_sample = c(1, 0.25, 1),
    data = gbti,
    quantiles = c(0.1, 0.9)
  )

  phase1_fit <- lme4::lmer(
    Response ~ Month + (Month | Patient),
    data = gbti,
    REML = TRUE
  )
  phase1_ranef <- lme4::ranef(phase1_fit)$Patient

  subject_id <- intersect(colnames(design$z_i), rownames(phase1_ranef))
  expect_equal(
    unname(design$z_i["intercept", subject_id]),
    unname(phase1_ranef[subject_id, "(Intercept)"]),
    tolerance = 1e-8
  )
  expect_equal(
    unname(design$z_i["slope", subject_id]),
    unname(phase1_ranef[subject_id, "Month"]),
    tolerance = 1e-8
  )

  vc_phase1 <- as.data.frame(lme4::VarCorr(phase1_fit))[,"sdcor"]
  expected_ests <- c(
    lme4::fixef(phase1_fit),
    log(as.numeric(vc_phase1[1:2])),
    log((1 + as.numeric(vc_phase1[3])) / (1 - as.numeric(vc_phase1[3]))),
    log(as.numeric(vc_phase1[4]))
  )

  expect_equal(design$xcol.phase1, c("(Intercept)", "Month"))
  expect_equal(unname(design$ests.phase1), unname(expected_ests), tolerance = 1e-8)
})

expect_bds_gbti_reference <- function(method, p_sample, quantiles, expected,
                                      coef_tol = 1e-6, vcov_tol = 1e-6)
{
  data(gbti)

  set.seed(202405)
  design <- bds(
    Response ~ Month | Patient,
    method,
    p_sample = p_sample,
    data = gbti,
    quantiles = quantiles
  )

  if (expected$code == 3L) {
    expect_warning(
      fit <- acml(Response ~ Month * Genotype, design),
      regexp = "last global step failed"
    )
  } else {
    fit <- acml(Response ~ Month * Genotype, design)
  }

  input <- bds_test_acml_input(
    design,
    as.character(unique(gbti$Patient[!is.na(gbti$Genotype)]))
  )
  w_function <- if (method == "bivariate") {
    "blup.bivariate"
  } else {
    paste0("blup.", method)
  }

  current <- LogLikeCAndScore2(
    expected$coefficients,
    y = input$y,
    x = input$x,
    z = input$z,
    id = input$id,
    w.function = rep(w_function, length(input$y)),
    cutpoints = input$cutpoints,
    SampProb = input$SampProb,
    xcol.phase1 = design$xcol.phase1,
    ests.phase1 = design$ests.phase1
  )

  expect_true(inherits(fit, "acml"))
  expect_equal(fit$Code, expected$code)
  expect_equal(length(unique(input$id)), expected$n_subjects)
  expect_equal(length(input$y), expected$n_rows)
  expect_equal(unname(design$ests.phase1), expected$phase1, tolerance = 1e-8)
  expect_equal(as.vector(design$cutpoints), expected$cutpoints, tolerance = 1e-8)
  expect_equal(coef(fit), expected$coefficients, tolerance = coef_tol)
  expect_equal(as.numeric(logLik(fit)), expected$logLik, tolerance = coef_tol)
  expect_equal(as.numeric(current), -expected$logLik, tolerance = coef_tol)
  expect_equal(attr(current, "gradient"), expected$gradient, tolerance = 1e-5)
  expect_equal(diag(fit$covariance), expected$covariance_diag, tolerance = vcov_tol)
  expect_equal(diag(fit$robcov), expected$robcov_diag, tolerance = vcov_tol)
}

test_that("BDS intercept estimates match fixed ACML.R gbti benchmark values",
{
  expect_bds_gbti_reference(
    "intercept",
    p_sample = c(1, 0.25, 1),
    quantiles = c(0.1, 0.9),
    expected = expected_bds_gbti_reference$intercept
  )
})

test_that("BDS slope estimates match fixed ACML.R gbti benchmark values",
{
  expect_bds_gbti_reference(
    "slope",
    p_sample = c(1, 0.25, 1),
    quantiles = c(0.1, 0.9),
    expected = expected_bds_gbti_reference$slope
  )
})

test_that("BDS bivariate estimates match fixed ACML.R gbti benchmark values",
{
  expect_bds_gbti_reference(
    "bivariate",
    p_sample = c(0.25, 1),
    quantiles = 0.8,
    expected = expected_bds_gbti_reference$bivariate,
    coef_tol = 1e-4,
    vcov_tol = 1e-4
  )
})

test_that("BDS slope likelihood matches the built-in ACML validation helper",
{
  data(gbti)

  set.seed(202405)
  design <- bds(
    Response ~ Month | Patient,
    "slope",
    p_sample = c(1, 0.25, 1),
    data = gbti,
    quantiles = c(0.1, 0.9)
  )

  input <- bds_test_acml_input(design, bds_test_reference_ids(gbti))
  params <- bds_test_params(input$x)
  w_function <- rep("blup.slope", length(input$y))

  current <- LogLikeCAndScore2(
    params,
    y = input$y,
    x = input$x,
    z = input$z,
    id = input$id,
    w.function = w_function,
    cutpoints = input$cutpoints,
    SampProb = input$SampProb,
    xcol.phase1 = design$xcol.phase1,
    ests.phase1 = design$ests.phase1
  )
  reference <- av_bds_LogLikeCAndScore2(
    params,
    y = input$y,
    x = input$x,
    z = input$z,
    id = input$id,
    w.function = w_function,
    cutpoints = input$cutpoints,
    SampProb = input$SampProb,
    Weights = rep(1, length(input$y)),
    xcol.phase1 = design$xcol.phase1,
    ests.phase1 = design$ests.phase1
  )

  expect_equal(as.numeric(current), as.numeric(reference), tolerance = 1e-8)
  expect_equal(
    unname(attr(current, "gradient")),
    unname(attr(reference, "gradient")),
    tolerance = 1e-8
  )
})

test_that("acml estimates from a BDS slope design match the built-in ACML validation helper",
{
  data(gbti)

  set.seed(202405)
  fit_data <- bds_test_fit_data(gbti)
  design <- bds(
    Response ~ Month | Patient,
    "slope",
    p_sample = c(1, 0.25, 1),
    data = fit_data,
    quantiles = c(0.1, 0.9)
  )

  input <- bds_test_acml_input(design, as.character(unique(fit_data$Patient)))
  init <- bds_test_params(input$x)

  current <- acml(
    Response ~ Month * Genotype,
    design,
    init = init
  )

  reference_data <- stats::model.frame(
    Response ~ Month * Genotype + Patient + wt_ref,
    data = design$data
  )

  reference <- av_bds_acml.lmem2(
    formula.fixed = Response ~ Month * Genotype,
    formula.random = ~ Month,
    data = reference_data,
    id = Patient,
    w.function = rep("blup.slope", nrow(reference_data)),
    InitVals = init,
    cutpoints = input$cutpoints,
    SampProb = input$SampProb,
    Weights = wt_ref,
    xcol.phase1 = design$xcol.phase1,
    ests.phase1 = design$ests.phase1
  )

  expect_equal(current$coefficients, reference$coefficients, tolerance = 1e-8)
  expect_equal(current$logLik, reference$logLik, tolerance = 1e-8)
  expect_equal(current$covariance, reference$covariance, tolerance = 1e-7)
  expect_equal(current$robcov, reference$robcov, tolerance = 1e-7)
})

test_that("BDS bivariate likelihood uses the bivariate BLUP correction",
{
  data(gbti)
  set.seed(202405)

  design <- bds(
    Response ~ Month | Patient,
    "bivariate",
    p_sample = c(0.25, 1),
    data = gbti,
    quantiles = 0.8
  )

  input <- bds_test_acml_input(design, bds_test_reference_ids(gbti, n_each = 2))
  params <- bds_test_params(input$x)

  expect_error(
    result <- LogLikeCAndScore2(
      params,
      y = input$y,
      x = input$x,
      z = input$z,
      id = input$id,
      w.function = rep("blup.bivariate", length(input$y)),
      cutpoints = input$cutpoints,
      SampProb = input$SampProb,
      xcol.phase1 = design$xcol.phase1,
      ests.phase1 = design$ests.phase1
    ),
    NA
  )
  expect_true(is.finite(as.numeric(result)))
  expect_true(all(is.finite(attr(result, "gradient"))))
})
