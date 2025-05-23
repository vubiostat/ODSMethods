library(mvtnorm)
library(ODSMethods)

generate_gbti <- function(N=1000, M=5, seed=1, w_function='bivariate')
{
  set.seed(seed)

  ancestry <- sample(c("Non-Hispanic", "Hispanic"), replace=TRUE, N)
  genotype <- rbinom(N, 1, 0.1+0.15*(ancestry == 'Hispanic'))

  bi <- rmvnorm(N, mean=c(0,0), sigma=
      matrix(c(4, 0.1, 0.1, 0.25), nrow=2, byrow=TRUE))

  gbti <- data.frame(
    Patient  = rep(1:N, each=M),
    Month    = rep(1:M, N),
    Ancestry = factor(rep(ancestry, each=M)),
    Response = NA,
    pSample  = NA,
    Genotype = rep(genotype==1, each=M)
  )

  gbti$Response <- with(gbti,
    10+0.5*Month-0.5*Genotype+0.25*Month*Genotype+0.5*(Ancestry=='Hispanic')+bi[Patient,1]+bi[Patient,2]*Month+rnorm(N*M))

  alpha <- 0.5*(1-sqrt(0.8))
  quantile <- if(w_function=='bivariate') c(alpha, 1-alpha) else c(0.2, 0.8)
  design <- ods(Response ~ Month|Patient,
                w_function,
                p_sample=c(1, 0.25, 1),
                data=gbti,
                quantiles=quantile)
  gbti$pSample <- rep(design$p_sample_i, each=M)

  # Sample it!
  gbti$Genotype[!rep(rbinom(N, 1, design$p_sample_i),each=M) == 1] <- NA

  gbti
}

gbti <- generate_gbti()
