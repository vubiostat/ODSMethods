context("ACML Continuous Likelihood Methods")

# This Tests the Package code versus the validated reference (for what we can test)
# At present this is limited to the logLik value.
# FIXME: Add Gradient Tests

data(gbti)

test_that("ACML continuous response mean",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'mean',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                gbti,
                init=rep(1, 8))
  )

  # Did it converge and return an object of the right type?
  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   1) # <- This will change based on optimization utilized

  valid <- validated_ll(gbti, coef(est), 'mean', c(0.1, 0.9), c(1, 0.25, 1))

  expect_close(logLik(est), valid)
})

test_that("ACML continuous response intercept",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'intercept',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                gbti,
                init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,  1) # <- This will change based on optimization utilized

  valid <- validated_ll(gbti, coef(est), 'intercept', c(0.1, 0.9), c(1, 0.25, 1))

  expect_close(logLik(est), valid)
})

test_that("ACML continuous response slope",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'slope',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
               design,
               gbti,
               init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,1) # This changes based on method
  valid <- validated_ll(gbti, coef(est), 'slope', c(0.1, 0.9), c(1, 0.25, 1))

  expect_close(logLik(est), valid)

})

test_that("ACML continuous response bivariate",
{
  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'bivariate',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                gbti,
                init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   1) # This changes based on method

  valid <- validated_ll(gbti, coef(est), 'bivar', c(0.1, 0.9), c(1, 0.25, 1))

  expect_close(logLik(est), valid)
})
