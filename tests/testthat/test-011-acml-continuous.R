context("ACML Continuous Likelihood Methods")

data(GroupByTimeInteraction)

data <- GroupByTimeInteraction

test_that("ACML continuous response intercept",
{
  cutpoints <- quantile(
    sapply(
      split(data, data$patient),
      function(d) coef(lm(response~month, d))[1]
    ),
    probs=c(0.1, 0.9)
  )

  design <- ods(response ~ (month|patient), data, c(1, 0.25, 1), cutpoints, 'intercept')

  est    <- acml(response ~ month*genotype, design)

  expect_equal(est$Code,   2)

  expect_close(est$Ests,   estimates_acml_intercept)

  expect_close(est$LogL,   logl_acml_intercept)

  expect_close(est$covar,  cv_acml_intercept)

  expect_close(est$robcov, rcv_acml_intercept)
})
