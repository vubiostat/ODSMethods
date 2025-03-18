context("ACML Continuous Design / Analysis")

# Comment out to validate values versus original published acml code
# NOTE: These are messed up and mis-specified, and only for point verification
#skip(message = "validated original acml code")

data(gbti)

test_that("acml handles character id",
{
  data <- gbti[gbti$Patient <= 200,]
  data$Patient <- as.character(data$Patient)

  expect_silent(
    design <- ods(Response ~ Month|Patient,
                  'mean',
                  p_sample=c(1, 0.25, 1),
                  data=data,
                  quantiles=c(0.1, 0.9))
  )

  expect_silent(
    est <- acml(Response ~ Month*Genotype,
                design,
                data,
                init=rep(1, 8))
  )

  expect_true(inherits(est, "acml"))
  expect_equal(est$Code,   2) # This changes based on method
  expect_close(coef(est),   estimates_acml_mean)
})

test_that("acml checks for ods mismatch",
{
  local_reproducible_output(width = 200)

  expect_silent(
    design <- ods(Response ~ Month|Patient,
                'mean',
                p_sample=c(1, 0.25, 1),
                data=gbti,
                subset=Patient <= 200,
                quantiles=c(0.1, 0.9))
  )

  expect_error(
    acml( Response ~ Month*Genotype,
          design,
          gbti,
          subset=Patient <= 300,
          init=rep(1, 8)),
    'Group variables provided to acml that were not part of design'
  )

  expect_error(
    acml( Response ~ Month*Genotype,
          design,
          gbti,
          subset=Patient >= 700,
          init=rep(1, 8)),
    'Group variables provided to acml that were not part of design'
  )
})
