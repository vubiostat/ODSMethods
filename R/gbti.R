#' Simulated Group by Time Interaction (gbti) dataset
#'
#' @name gbti
#' @docType data
#' @description
#' This is a simulated data set with known ``truth''. The dataset consists of
#' 402 patients and a measured response to a hypothetical drug over the course
#' of 5 months. Further, some ethnicity information is available if
#' Hispanic origin or not. An expensive genotype analysis exists, and to
#' minimize exposure a biased sample is performed. Patients that are in the
#' outer 20% of responses for intercept and slope are sampled for the genotyping
#' with a probability of 1. Patients that are in the inner 80% of responses for
#' intercept and slope are sampled with a probability of 0.25 for genotyping.
#' Thus on average 304 patients of the 1000 would be genotyped.
#'
#' @details
#' The model was generated using the following:
#' \deqn{Y_{ij} = 10 + 0.5 \cdot Month - 0.5 \cdot Genotype  + 0.25 \cdot Month
#' \cdot Genotype + 0.5 \cdot Ancestry + b_{0i} + b_{1i} \cdot Month + \epsilon_{ij}}
#'
#' where
#'
#' \deqn{b_{0i} \sim \mathcal{N}(\mu=0, \sigma^2_0=4)}
#' \deqn{b_{1i} \sim \mathcal{N}(\mu=0, \sigma^2_1=0.25)}
#' \deqn{\rho(b_{0i}, b_{1i}) = 0.1}
#' \deqn{\epsilon_{ij} \sim \mathcal{N}(\mu=0, \sigma^2=1)}
#'
#' Further,
#'
#' \deqn{\pi(genotype_i = 1 | ancestry_i) = 0.1 + 0.15 \cdot I(ancestry_i=1)}
#'
#' Assume that funding to genotype only 400 individuals is
#' available, if 200 in the outer bivariate 20\% are genotyped at 100%, this
#' means that the remaining 800 can only be genotyped at 25\%. The
#' quantiles of 5.28\% and 94.72\% give the desired cut points,
#' \eqn{\frac{1}{2}(1-\sqrt{0.80})}.
#'
#' @usage data(gbti)
#' @return Patient responses, with stratified sampling of genotype in a biased manner based on outcome.
#' @examples
#'   data(gbti)
NULL


