context("ACML Validated Continuous Likelihood Methods")

data(gbti)

data <- gbti

b_i <-sapply(
        split(data, data$Patient),
        function(d) c(mean(d$Response),coef(lm(Response~Month, d)))
      )
rownames(b_i) <- c("mean", "intercept", "slope")

data <- data[!is.na(data$Genotype),]  # Drop NA data for these tests

# Create the design matrices
d_mat <- cbind(rep(1, nrow(data)),
               data[,c("Month", "Genotype")])
d_mat$g_t <- d_mat$Month*d_mat$Genotype

y  <- matrix(data$Response,ncol = 1)
x  <- as.matrix(d_mat)
z  <- as.matrix(cbind(rep(1, length(data$Month)), data$Month))
id <- matrix(data$Patient,ncol = 1)

test_that("Validation ACML 'mean' reference works",
{
  cuts <- quantile(b_i['mean',], probs=c(0.1, 0.9))

  result <- av_LogLikeAndScore(
    estimates_acml_mean, y, x, z, id, 'mean', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_mean)
  expect_close(attr(result, "gradient"), gradient_acml_mean)
})

test_that("Validation ACML 'intercept' reference works",
{
  cuts <- quantile(b_i['intercept',], probs=c(0.1, 0.9))

  result <- av_LogLikeAndScore(
    estimates_acml_intercept, y, x, z, id, 'intercept', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_intercept)
  expect_close(attr(result, "gradient"), gradient_acml_intercept)
})

test_that("Validation ACML 'slope' reference works",
{
  cuts <- quantile(b_i['slope',], probs=c(0.1, 0.9))

  result <- av_LogLikeAndScore(
    estimates_acml_slope, y, x, z, id, 'slope', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_slope)
  expect_close(attr(result, "gradient"), gradient_acml_slope)
})

test_that("Validation ACML 'bivariate' reference works",
{
  cuts <- c(
    quantile(b_i['intercept',], probs=c(0.1, 0.9)),
    quantile(b_i['slope',], probs=c(0.1, 0.9))
  )

  result <- av_LogLikeAndScore(
    estimates_acml_bivar, y, x, z, id, 'bivar', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_bivar)
  expect_close(attr(result, "gradient"), gradient_acml_bivar)
})
