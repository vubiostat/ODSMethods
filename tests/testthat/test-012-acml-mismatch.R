context("ACML Continuous Design / Analysis")

# Comment out to validate values versus original published acml code
# NOTE: These are messed up and mis-specified, and only for point verification
#skip(message = "validated original acml code")

data(gbti)

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
