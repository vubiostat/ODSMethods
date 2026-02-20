
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
    p0   <- cuts[1]
    ints <- b_i["intercept", ]
    slps <- b_i["slope",     ]

    q1  <- 0.99
    Del <- 1

    ## marginal joint inside around p0
    while (Del > 0.003 && q1 > 0.5) {
      q1 <- q1 - 0.001

      I_low  <- as.numeric(stats::quantile(ints, probs = 1 - q1))
      I_high <- as.numeric(stats::quantile(ints, probs = q1))
      S_low  <- as.numeric(stats::quantile(slps, probs = 1 - q1))
      S_high <- as.numeric(stats::quantile(slps, probs = q1))

      inside <- (ints > I_low & ints < I_high &
                   slps > S_low & slps < S_high)

      Del <- abs(mean(inside) - p0)
    }

    c(I_low, I_high, S_low,  S_high)

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
