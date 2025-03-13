context("ACML Validated Continuous Likelihood Methods")

# Comment out to validate values versus original published acml code
# NOTE: These are messed up and mis-specified, and only for point verification
#skip(message = "validated original acml code")

data(gbti)

data <- gbti[gbti$Patient <= 200,] # Only use patients 1:200

b_i <-
    sapply(
      split(data, data$Patient),
      function(d) c(mean(d$Response),coef(lm(Response~Month, d)))
    )
rownames(b_i) <- c("mean", "intercept", "slope")


data <- data[!is.na(data$Genotype),]  # Drop NA data for these tests

# Create a design matrix
d_mat <- cbind(rep(1, nrow(data)),
               data[,c("Month", "Genotype")])
d_mat$g_t <- d_mat$Month*d_mat$Genotype

y  <- matrix(data$Response,ncol = 1)
x  <- as.matrix(d_mat)
z  <- as.matrix(cbind(rep(1, length(data$Month)), data$Month))
id <- matrix(data$Patient,ncol = 1)

test_that("Validation ACML 'mean' reference works",
{
  # Note: This is an incorrect practice as it's on the post-sampled data
  # However, for a point test of the method it's fine.
  cuts <- quantile(b_i['mean',], probs=c(0.1, 0.9))

  validated_acml_mean <-
    acml_validated(
      y=y, x=x, z=z, id=id,
      w.function = "mean",
      InitVals = c(1,1,1,1,1,1,1,1),
      cutpoints = cuts,
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
  cuts <- quantile(b_i['intercept',], probs=c(0.1, 0.9))

  validated_acml_intercept <-
    acml_validated(
      y=y, x=x, z=z, id=id,
      w.function = "intercept",
      InitVals = c(1,1,1,1,1,1,1,1),
      cutpoints = cuts,
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
  cuts <- quantile(b_i['slope',], probs=c(0.1, 0.9))

  validated_acml_slope <-
    acml_validated(
      y=y, x=x, z=z, id=id,
      w.function = "slope",
      InitVals = c(1,1,1,1,1,1,1,1),
      cutpoints = cuts,
      SampProb = c(1, 0.25, 1)
    )

  expect_equal(validated_acml_slope$Code,   2)
  expect_close(validated_acml_slope$Ests,   estimates_acml_slope)
  expect_close(validated_acml_slope$LogL,   logl_acml_slope)
  expect_close(validated_acml_slope$covar,  cv_acml_slope)
  expect_close(validated_acml_slope$robcov, rcv_acml_slope)
})

test_that("Validation ACML 'bivariate' reference works",
{
  cuts <- c(
    quantile(b_i['intercept',], probs=c(0.1, 0.9)),
    quantile(b_i['slope',], probs=c(0.1, 0.9))
  )

  validated_acml_bivar <-
    acml_validated(
      y=y, x=x, z=z, id=id,
      w.function = "bivar",
      InitVals = c(1,1,1,1,1,1,1,1),
      cutpoints = cuts,
      SampProb = c(1, 0.25, 1)
    )

  expect_equal(validated_acml_bivar$Code,   2)
  expect_close(validated_acml_bivar$Ests,   estimates_acml_bivar)
  expect_close(validated_acml_bivar$LogL,   logl_acml_bivar)
  expect_close(validated_acml_bivar$covar,  cv_acml_bivar)
  expect_close(validated_acml_bivar$robcov, rcv_acml_bivar)
})
