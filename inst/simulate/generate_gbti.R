library(mvtnorm)
library(ODSMethods)

set.seed(1)

ancestry <- sample(c("Non-Hispanic", "Hispanic"), replace=TRUE, 1000)
genotype <- rbinom(1000, 1, 0.1+0.15*(ancestry == 'Hispanic'))

bi <- rmvnorm(1000, mean=c(0,0), sigma=
    matrix(c(4, 0.1, 0.1, 0.25), nrow=2, byrow=TRUE))

gbti <- data.frame(
  Patient  = rep(1:1000, each=5),
  Month    = rep(1:5, 5),
  Ancestry = factor(rep(ancestry, each=5)),
  Response = NA,
  pSample  = NA,
  Genotype = rep(genotype==1, each=5)
)

gbti$Response <- with(gbti,
  10+0.5*Month-0.5*Genotype+0.25*Month*Genotype+0.5*(Ancestry=='Hispanic')+bi[Patient,1]+bi[Patient,2]*Month+rnorm(5000))

alpha <- 0.5*(1-sqrt(0.8))
design <- ods(Response ~ Month|Patient,
                  'bivariate',
                  p_sample=c(1, 0.25, 1),
                  data=gbti,
                  quantiles=c(alpha, 1-alpha))
gbti$pSample <- rep(design$p_sample_i, each=5)

gbti$Genotype[!rep(rbinom(1000, 1, design$p_sample_i),each=5) == 1] <- NA

