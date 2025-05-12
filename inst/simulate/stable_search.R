source("inst/simulate/generate_gbti.R")

search <- lapply(1:10,function(x)
{
  data <- generate_gbti(N=250, seed=x)

  alpha <- 0.5*(1-sqrt(0.8))

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

  ODSMethods:::acml_hessian(c(10, 0.5, -0.5, 0.25, log(2), log(0.5), atanh(0.1), log(1)),
    y, x, z, id,c(alpha, 1-alpha), 'bivariate')
})
