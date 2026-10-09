context("MI Validated Type-1 Predictive Mean Matching")

expect_validated_mi <- function(method, expected = NULL, data = validated_mi_data())
{
  design <- validated_mi_design(method, data)
  reference <- validated_mi(design)
  dat <- design$data
  dat$grp[dat$sampled == 0] <- NA
  fit <- MI(y ~ grp * time + conf, design,
            exposure = grp, exposure_model = ~ conf, data = dat,
            M = 3, seed = 37, keep_completed = TRUE)

  expect_true(inherits(fit, "MI"))
  expect_true(all(reference$codes %in% c(1, 2)))
  expect_close(coef(fit), reference$coefficients, tol = 1e-6)
  expect_close(vcov(fit), reference$covariance, tol = 1e-6)
  expect_close(fit$coef_mi, reference$coef_mi, tol = 1e-6)
  expect_close(fit$U_bar, reference$U_bar, tol = 1e-6)
  expect_close(fit$B, reference$B, tol = 1e-6)
  for (m in 1:3) {
    expect_equal(fit$completed_data[[m]]$grp, reference$completed_data[[m]]$grp)
    expect_close(fit$cov_mi[[m]], reference$cov_mi[[m]], tol = 1e-6)
  }
  if (!is.null(expected)) {
    expect_close(coef(fit), expected$coefficients, tol = 1e-6)
    expect_close(diag(vcov(fit)), expected$variance, tol = 1e-6)
  }
}

test_that("MI continuous response mean",
{
  expect_validated_mi("mean", expected_mi_reference$mean)
})

test_that("MI continuous response intercept",
{
  expect_validated_mi("intercept", expected_mi_reference$intercept)
})

test_that("MI continuous response slope",
{
  expect_validated_mi("slope", expected_mi_reference$slope)
})

test_that("MI continuous response bivariate",
{
  expect_validated_mi("bivariate", expected_mi_reference$bivariate)
})

test_that("MI supports continuous, binary, and ordinal exposures",
{
  for (exposure_type in c("continuous", "binary", "ordinal")) {
    dat <- validated_mi_data(continuous = exposure_type == "continuous",
                              exposure_type = exposure_type)
    expect_validated_mi("slope", data = dat)
  }
})
