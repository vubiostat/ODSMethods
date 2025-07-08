context("ACML Validated Continuous Likelihood Methods")


test_that("Validation ACML 'mean' reference works",
{
  result <- validated_ll(data, estimates_acml_mean, 'mean')

  expect_close(result,                   logl_acml_mean)
  expect_close(attr(result, "gradient"), gradient_acml_mean)
})

test_that("Validation ACML 'intercept' reference works",
{
  result <- validated_ll(data, estimates_acml_intercept, 'intercept')

  expect_close(result,                   logl_acml_intercept)
  expect_close(attr(result, "gradient"), gradient_acml_intercept)
})

test_that("Validation ACML 'slope' reference works",
{
  result <- validated_ll(data, estimates_acml_slope, 'slope')

  expect_close(result,                   logl_acml_slope)
  expect_close(attr(result, "gradient"), gradient_acml_slope)
})

test_that("Validation ACML 'bivariate' reference works",
{
  result <- validated_ll(data, estimates_acml_bivar, 'bivar')

  expect_close(result,                   logl_acml_bivar)
  expect_close(attr(result, "gradient"), gradient_acml_bivar)
})
