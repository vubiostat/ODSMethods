context("ACML Continuous Likelihood Methods")

data(gbti)


test_that("ACML continuous response mean",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'mean',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  subset=Patient <= 200,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                gbti,
                subset=Patient <= 200,
                init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   2) # This changes based on method
  expect_close(coef(est),   estimates_acml_mean)
  expect_close(logLik(est), logl_acml_mean)
  expect_close(vcov(est),   cv_acml_mean)
  expect_close(robcov(est), rcv_acml_mean)

})

test_that("ACML continuous response intercept",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'intercept',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  subset=Patient <= 200,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                gbti,
                subset=Patient <= 200,
                init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   2) # This changes based on method
  expect_close(coef(est),   estimates_acml_intercept)
  expect_close(logLik(est), logl_acml_intercept)
  expect_close(vcov(est),   cv_acml_intercept)
  expect_close(robcov(est), rcv_acml_intercept)
})

test_that("ACML continuous response slope",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'slope',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  subset=Patient <= 200,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                gbti,
                subset=Patient <= 200,
                init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   2) # This changes based on method
  expect_close(coef(est),   estimates_acml_slope)
  expect_close(logLik(est), logl_acml_slope)
  expect_close(vcov(est),   cv_acml_slope)
  expect_close(robcov(est), rcv_acml_slope)

})

test_that("ACML continuous response bivariate",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'bivariate',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  subset=Patient <= 200,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                gbti,
                subset=Patient <= 200,
                init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   2) # This changes based on method
  expect_close(coef(est),   estimates_acml_bivar)
  expect_close(logLik(est), logl_acml_bivar)
  expect_close(vcov(est),   cv_acml_bivar)
  expect_close(robcov(est), rcv_acml_bivar)
})
