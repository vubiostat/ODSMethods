context("WL Validated Continuous Likelihood Methods")

# This Tests the Package code versus the validated reference (for what we can test)
# At present this is limited to the logLik value.

data(gbti)

validated_wl_input <- function(design, formula)
{
  all_vars_formula <- as.formula(
    paste(
      design$response, "~",
      paste(unique(c(all.vars(formula),
                     design$time,
                     design$id,
                     design$weights)),
            collapse = " + ")
    )
  )

  mf <- model.frame(
    all_vars_formula,
    data = design$data[design$data$sampled == 1, , drop = FALSE],
    drop.unused.levels = TRUE
  )

  y <- matrix(mf[[design$response]], ncol = 1)
  fixed.mf <- model.frame(formula, mf)
  x <- model.matrix(formula, fixed.mf)
  z <- as.matrix(data.frame(
    `(Intercept)` = 1,
    mf[, design$time, drop = FALSE],
    check.names = FALSE
  ))
  id <- matrix(mf[[design$id]], ncol = 1)
  weights <- mf[[design$weights]]

  list(y = y, x = x, z = z, id = id, weights = weights)
}

validated_wl <- function(fit, design, formula)
{
  input <- validated_wl_input(design, formula)

  av_LogLikeAndScore(
    coef(fit),
    input$y,
    input$x,
    input$z,
    input$id,
    "mean",
    c(0, 1),
    c(1, 1, 1),
    1 / input$weights,
    NA
  )
}

package_wl <- function(fit, design, formula)
{
  input <- validated_wl_input(design, formula)

  LogLikeCAndScoreWL(
    coef(fit),
    y = input$y,
    x = input$x,
    z = input$z,
    id = input$id,
    Weights = input$weights
  )
}

expect_validated_wl <- function(method, p_sample, quantiles, expected)
{
  set.seed(20260813)
  design <- ods(Response ~ Month|Patient,
                method,
                p_sample = p_sample,
                data = gbti,
                quantiles = quantiles,
                weights = "pSample")

  fit <- WL(Response ~ Month * Genotype, design)

  expect_true(inherits(fit, "WL"))
  expect_equal(fit$Code, 1)
  expect_close(coef(fit), expected$coefficients, tol = 1e-6)

  package_result <- package_wl(fit, design, Response ~ Month * Genotype)
  reference_result <- validated_wl(fit, design, Response ~ Month * Genotype)

  expect_close(logLik(fit), expected$logLik, tol = 1e-6)
  expect_close(attr(package_result, "gradient"), expected$gradient, tol = 1e-6)
  expect_close(package_result, reference_result)
  expect_close(attr(package_result, "gradient"),
               attr(reference_result, "gradient"))
  expect_close(logLik(fit), -reference_result)
}

test_that("WL continuous response mean",
{
  expect_validated_wl("mean", c(1, 0.25, 1), c(0.1, 0.9),
                      expected_wl_reference$mean)
})

test_that("WL continuous response intercept",
{
  expect_validated_wl("intercept", c(1, 0.25, 1), c(0.1, 0.9),
                      expected_wl_reference$intercept)
})

test_that("WL continuous response slope",
{
  expect_validated_wl("slope", c(1, 0.25, 1), c(0.1, 0.9),
                      expected_wl_reference$slope)
})

test_that("WL continuous response bivariate",
{
  expect_validated_wl("bivariate", c(1, 0.25), 0.1,
                      expected_wl_reference$bivariate)
})

test_that("WL agrees with lme4 under complete sampling",
{
  full_data <- gbti[!is.na(gbti$Genotype), ]
  full_data$wt <- 1

  set.seed(20260821)
  design <- ods(Response ~ Month|Patient,
                "mean",
                p_sample = c(1, 1, 1),
                data = full_data,
                quantiles = c(0.1, 0.9),
                weights = "wt")

  fit_wl <- WL(Response ~ Month * Genotype, design)

  lme4_data <- design$data[design$data$sampled == 1, , drop = FALSE]
  fit_lme4 <- lme4::lmer(Response ~ Month * Genotype + (Month|Patient),
                         data = lme4_data,
                         weights = wt,
                         REML = FALSE)

  lme4_re <- as.data.frame(lme4::VarCorr(fit_lme4))$sdcor

  expect_close(fixef(fit_wl), lme4::getME(fit_lme4, "beta"), tol = 1e-6)
  expect_close(ranef(fit_wl, transform = TRUE), lme4_re, tol = 1e-4)
  expect_close(logLik(fit_wl), stats::logLik(fit_lme4), tol = 1e-6)
})
