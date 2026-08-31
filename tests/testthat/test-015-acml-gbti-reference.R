context("ACML gbti benchmark values")

data(gbti)

acml_gbti_input <- function(fit, design, formula)
{
  ftxt <- paste(deparse(design$formula), collapse = "")
  parts <- strsplit(ftxt, "\\|", perl = TRUE)[[1]]
  f_left <- stats::as.formula(trimws(parts[1]))
  x_vars <- all.vars(stats::update(f_left, 0 ~ .))

  all_vars <- unique(c(all.vars(formula), x_vars, design$id))
  all_vars_formula <- as.formula(
    paste(design$response, "~",
          paste(all_vars[all_vars != design$response], collapse = " + "))
  )

  mf <- model.frame(all_vars_formula, data = design$data,
                    drop.unused.levels = TRUE)
  fixed.mf <- model.frame(formula, mf)
  args <- attr(fit, "args")

  list(
    y = unlist(mf[, design$response, drop = FALSE]),
    x = model.matrix(formula, fixed.mf),
    z = as.matrix(data.frame(
      `(Intercept)` = 1,
      mf[, x_vars, drop = FALSE],
      check.names = FALSE
    )),
    id = args$id,
    w.function = args$w.function,
    cutpoints = args$cutpoints,
    SampProb = args$SampProb
  )
}

acml_gbti_gradient <- function(fit, design, formula)
{
  input <- acml_gbti_input(fit, design, formula)

  subjectData <- CreateSubjectData(
    id = input$id,
    y = input$y,
    x = input$x,
    z = input$z,
    SampProb = input$SampProb,
    cutpoints = input$cutpoints,
    w.function = input$w.function,
    xcol.phase1 = design$xcol.phase1,
    ests.phase1 = design$ests.phase1
  )

  out <- LogLikeCAndScore2(
    coef(fit),
    y = input$y,
    x = input$x,
    z = input$z,
    id = input$id,
    w.function = input$w.function,
    cutpoints = input$cutpoints,
    SampProb = input$SampProb,
    xcol.phase1 = design$xcol.phase1,
    ests.phase1 = design$ests.phase1,
    subjectData = subjectData
  )

  attr(out, "gradient")
}

expect_acml_gbti_reference <- function(method, p_sample, quantiles, expected)
{
  set.seed(20260831)
  design <- ods(Response ~ Month|Patient,
                method,
                p_sample = p_sample,
                data = gbti,
                quantiles = quantiles)

  fit <- acml(Response ~ Month * Genotype, design)

  expect_true(inherits(fit, "acml"))
  expect_equal(fit$Code, 1)
  input <- acml_gbti_input(fit, design, Response ~ Month * Genotype)

  expect_equal(length(unique(input$id)), expected$n_subjects)
  expect_equal(length(input$y), expected$n_rows)
  expect_close(as.vector(design$cutpoints), expected$cutpoints, tol = 1e-8)
  expect_close(coef(fit), expected$coefficients, tol = 1e-6)
  expect_close(logLik(fit), expected$logLik, tol = 1e-6)
  expect_close(acml_gbti_gradient(fit, design, Response ~ Month * Genotype),
               expected$gradient, tol = 1e-5)
}

test_that("ACML gbti mean estimates match fixed ACML.R benchmark values",
{
  expect_acml_gbti_reference("mean", c(1, 0.25, 1), c(0.1, 0.9),
                             expected_acml_gbti_reference$mean)
})

test_that("ACML gbti intercept estimates match fixed ACML.R benchmark values",
{
  expect_acml_gbti_reference("intercept", c(1, 0.25, 1), c(0.1, 0.9),
                             expected_acml_gbti_reference$intercept)
})

test_that("ACML gbti slope estimates match fixed ACML.R benchmark values",
{
  expect_acml_gbti_reference("slope", c(1, 0.25, 1), c(0.1, 0.9),
                             expected_acml_gbti_reference$slope)
})

test_that("ACML gbti bivariate estimates match fixed ACML.R benchmark values",
{
  expect_acml_gbti_reference("bivariate", c(1, 0.25), 0.1,
                             expected_acml_gbti_reference$bivariate)
})
