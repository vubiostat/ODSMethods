
# This function takes data and creates the design objects in the validated
# direct manner and calls the validated reference.
# It assumes the reference data matches the columns of gbti

validated_ll <- function(data, est, method, cuts=c(0.1, 0.9), probs=c(1, 0.25, 1))
{
  b_i <-sapply(
    split(data, data$Patient),
    function(d) c(mean(d$Response),coef(lm(Response~Month, d)))
  )
  rownames(b_i) <- c("mean", "intercept", "slope")

  cuts <- if(method == 'bivar')
  {
    c(
      quantile(b_i['intercept',], probs=cuts),
      quantile(b_i['slope',],     probs=cuts)
    ) # FIXME: this needs to be changed. The bivariate design is not a simple cut of intercept and slope
  } else quantile(b_i[method,], probs=cuts)

  data <- data[!is.na(data$Genotype),]  # Drop NA data for these tests

  # Create the design matrices
  d_mat <- cbind(rep(1, nrow(data)),
                 data[,c("Month", "Genotype")])
  d_mat$g_t <- d_mat$Month*d_mat$Genotype

  y  <- matrix(data$Response,ncol = 1)
  x  <- as.matrix(d_mat)
  z  <- as.matrix(cbind(rep(1, length(data$Month)), data$Month))
  id <- matrix(data$Patient,ncol = 1)

  -1.0*av_LogLikeAndScore(
    est, y, x, z, id, method, cuts, c(1, 0.25, 1),
    rep(1, length(y)), NA)
}
