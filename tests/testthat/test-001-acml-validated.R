context("ACML Validated Continuous Likelihood Methods")

data(gbti)

data <- gbti
data <- data[!is.na(data$genotype),]  # Drop NA data for these tests

# Create a design matrix
d_mat <- cbind(rep(1, nrow(data)),
               data[,c("month", "genotype")])
d_mat$g_t <- d_mat$month*d_mat$genotype

y  <- matrix(data$response,ncol = 1)
x  <- as.matrix(d_mat)
z  <- as.matrix(cbind(rep(1, length(data$month)), data$month))
id <- matrix(data$patient,ncol = 1)

test_that("Validation ACML 'mean' reference works",
{
  means <- quantile(
    aggregate(data$response, list(data$patient), FUN=mean)$x,
    probs=c(0.1, 0.9)
  )

  validated_acml_mean <-
    acml_validated(
      y=y, x=x, z=z, id=id,
      w.function = "mean",
      InitVals = c(1,1,1,1,1,1,1,1),
      cutpoints = means,
      SampProb = c(1, 0.25, 1)
    )

  expect_equal(validated_acml_mean$Code,   2)
  expect_close(validated_acml_mean$Ests,   estimates_acml_mean)
  expect_close(validated_acml_mean$LogL,   logl_acml_mean)
  expect_close(validated_acml_mean$covar,  cv_acml_mean)
  expect_close(validated_acml_mean$robcov, rcv_acml_mean)
})

test_that("Validation ACML 'intercept' reference works",
{
  intercepts <- quantile(
    sapply(
      split(data, data$patient),
      function(d) coef(lm(response~month, d))[1]
    ),
    probs=c(0.1, 0.9)
  )

  validated_acml_intercept <-
    acml_validated(
      y=y, x=x, z=z, id=id,
      w.function = "intercept",
      InitVals = c(1,1,1,1,1,1,1,1),
      cutpoints = intercepts,
      SampProb = c(1, 0.25, 1)
    )

  expect_equal(validated_acml_intercept$Code,   2)
  expect_close(validated_acml_intercept$Ests,   estimates_acml_intercept)
  expect_close(validated_acml_intercept$LogL,   logl_acml_intercept)
  expect_close(validated_acml_intercept$covar,  cv_acml_intercept)
  expect_close(validated_acml_intercept$robcov, rcv_acml_intercept)
})

test_that("Validation ACML 'slope' reference works",
{
  slopes <- quantile(
    sapply(
      split(data, data$patient),
      function(d) coef(lm(response~month, d))[2]
    ),
    probs=c(0.1, 0.9)
  )

  validated_acml_slope <-
    acml_validated(
      y=y, x=x, z=z, id=id,
      w.function = "slope",
      InitVals = c(1,1,1,1,1,1,1,1),
      cutpoints = slopes,
      SampProb = c(1, 0.5, 1)
    )

  expect_equal(validated_acml_slope$Code,   2)
  expect_close(validated_acml_slope$Ests,   estimates_acml_slope)
  expect_close(validated_acml_slope$LogL,   logl_acml_slope)
  expect_close(validated_acml_slope$covar,  cv_acml_slope)
  expect_close(validated_acml_slope$robcov, rcv_acml_slope)
})
