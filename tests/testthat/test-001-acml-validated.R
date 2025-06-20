context("ACML Validated Continuous Likelihood Methods")

# Comment out to validate values versus original published acml code
# NOTE: These are messed up and mis-specified, and only for point verification
#skip(message = "validated original acml code")

data(gbti)

data <- gbti

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
  # Expected values from reference
  estimates_acml_mean <-
    c(10.2355307756529, 0.497701195031362, -0.869402426613556, 0.234975440487394,
      0.804140810380354, -0.646309664693433, -0.857060994416562, -0.0952586041055383
    )
  logl_acml_mean     <- 1113.47834702023
  gradient_acml_mean <-c(-6.40822068465052e-07,
                         -1.02993776822657e-05, 3.4465476994594e-07, -2.0906257773845e-06,
                         -2.38080078014845e-07, 5.92712023121566e-06, 2.1020533366353e-06,
                         5.32940419211524e-06)

  cuts <- quantile(b_i['mean',], probs=c(0.1, 0.9))

  result <- LogLikeAndScore(
    estimates_acml_mean, y, x, z, id, 'mean', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_mean)
  expect_close(attr(result, "gradient"), gradient_acml_mean)
})

test_that("Validation ACML 'intercept' reference works",
{
  estimates_acml_intercept <-
    c(10.2571120507431, 0.505280246496726, -0.535325022821004, 0.198390388964521,
      0.630422450934662, -0.640158610225613, -0.0830362372320438, -0.0952586088753488
    )

  logl_acml_intercept <- 1129.90154093494
  gradient_acml_intercept <- c(-2.28652900169379e-07,
                               -3.31806911280097e-06, -2.04208166487696e-07, -1.47732674804502e-06,
                               4.18834839483793e-08, -4.97512139853957e-07, -3.68459503710184e-08,
                               7.41434245909591e-07)

  cuts <- quantile(b_i['intercept',], probs=c(0.1, 0.9))

  result <- LogLikeAndScore(
    estimates_acml_intercept, y, x, z, id, 'intercept', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_intercept)
  expect_close(attr(result, "gradient"), gradient_acml_intercept)
})

test_that("Validation ACML 'slope' reference works",
{
  estimates_acml_slope <-
    c(10.4876464692288, 0.440935051119677, -1.80251674681208, 0.282185161953021,
      1.0203795, -0.787915936837711, 0.26544380669411, 0.0343661690629439
    )
  logl_acml_slope <- 1139.62730702754
  gradient_acml_slope <- c(0.518651679002783, -23.1486323458273,
                           -1.61843833442785, 11.546578348304, 17.9892990412853, 61.5588249249227,
                           11.1289097518247, 114.264504032877)
  cuts <- quantile(b_i['slope',], probs=c(0.1, 0.9))

  result <- LogLikeAndScore(
    estimates_acml_slope, y, x, z, id, 'slope', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_slope)
  expect_close(attr(result, "gradient"), gradient_acml_slope)
})


test_that("Validation ACML 'bivariate' reference works",
{
  estimates_acml_bivar <-
    c(10.7655571285408, 0.246129113385224, -2.98989684370386, 0.801707586626488,
      1.27691825783006, -0.137203695930768, -0.0408259494816768, 0.0343661614845943
    )
  logl_acml_bivar <- 1214.58723417167
  gradient_acml_bivar <- c(0.430732452414889, -17.8598196167201,
                           -0.98725870359207, 2.387990523425, 5.48468228857626, 17.3081935840506,
                           8.2789001784875, 92.7896770053413)

  cuts <- c(
    quantile(b_i['intercept',], probs=c(0.1, 0.9)),
    quantile(b_i['slope',], probs=c(0.1, 0.9))
  )

  result <- LogLikeAndScore(
    estimates_acml_bivar, y, x, z, id, 'bivar', cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)

  expect_close(result,                   logl_acml_bivar)
  expect_close(attr(result, "gradient"), gradient_acml_bivar)
})
