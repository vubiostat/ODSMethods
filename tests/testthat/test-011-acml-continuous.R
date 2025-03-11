context("ACML Continuous Likelihood Methods")

data(GroupByTimeInteraction)

data <- GroupByTimeInteraction

test_that("ACML continuous response intercept",
{
  design <- ods(response ~ month|patient,
                'intercept',
                p_sample=c(1, 0.25, 1),
                data=data,
                quantiles=c(0.1, 0.9))

  est    <- acml(response ~ month*genotype, design, data, init=rep(1, 8))

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   2) # This changes based on method
  expect_close(coef(est),   estimates_acml_intercept)
  expect_close(est$LogL,    logl_acml_intercept)
  expect_close(vcov(est),   rcv_acml_intercept)

  #expect_true("residuals" %in% names(est))
  #expect_true("fitted.values" %in% names(est))
  #expect_true("rank" %in% names(est))
  #expect_true("df.residual" %in% names(est))

})
