// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix vi_calc(
    const Rcpp::NumericMatrix& zi,
    const Rcpp::NumericVector& sigma_vc,
    const Rcpp::NumericVector& rho_vc,
    const Rcpp::NumericVector& sigma_e
)
{
  // nrowzi    <- nrow(zi)
  int nrowzi = zi.nrow();
  // ncolzi <- ncol(zi) ## make sure this equals length(sigma.vc)
  int ncolzi = zi.ncol();
  // nERRsd <- length(sigma.e)
  int nERRsd = sigma_e.size();

  // Question: This is a vector in R?
  // SDMat.RE  <- diag(sigma.vc)
  // Diagonal matrix of random-effect standard deviations
  Rcpp::NumericMatrix SDMat_RE(ncolzi, ncolzi);
  for (int i = 0; i < ncolzi; ++i)
    SDMat_RE(i, i) = sigma_vc[i];

  // Construct correlation matrix
  Rcpp::NumericMatrix CorMat_RE(ncolzi, ncolzi);
  // b         <- matrix(0,ncolzi,ncolzi)
  // b[lower.tri(b, diag=FALSE)] <- rho.vc
  // CorMat.RE <- t(b)+b+diag(rep(1,ncolzi))

  int rho_idx = 0;

  for (int i = 0; i < ncolzi; ++i)
  {
    CorMat_RE(i, i) = 1.0;

    for (int j = 0; j < i; ++j)
    {
      CorMat_RE(i, j) = rho_vc[rho_idx];
      CorMat_RE(j, i) = rho_vc[rho_idx];

      ++rho_idx;
    }
  }

  // Covariance matrix of random effects:
  //
  // CovMat.RE <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
  Rcpp::NumericMatrix CovMat_RE(ncolzi, ncolzi);

  for (int i = 0; i < ncolzi; ++i)
  {
    for (int j = 0; j < ncolzi; ++j)
    {
      CovMat_RE(i, j) = sigma_vc[i]*CorMat_RE(i, j)*sigma_vc[j];
    }
  }

  // V_i = Z_i D Z_i' + sigma_e^2 I
  Rcpp::NumericMatrix V(nrowzi, nrowzi);

  // R version: V_i = zi %*% CovMat.RE %*% t(zi)
  for (int i = 0; i < nrowzi; ++i)
  {
    for (int j = 0; j < nrowzi; ++j)
    {
      double value = 0.0;

      for (int k = 0; k < ncolzi; ++k)
      {
        for (int l = 0; l < ncolzi; ++l)
        {
          value += zi(i, k)*CovMat_RE(k, l)*zi(j, l);
        }
      }

      V(i, j) = value;
    }
  }

  // Add measurement-error variance to the diagonal.
  //
  // R version: zi %*% CovMat.RE %*% t(zi) + diag(rep(sigma.e^2, each=nrowzi/nERRsd))
  // diag(rep(sigma.e^2, each=nrowzi/nERRsd))
  int observations_per_sigma = nrowzi / nERRsd;

  for (int i = 0; i < nrowzi; ++i)
  {
    int sigma_idx = i / observations_per_sigma;

    V(i, i) += sigma_e[sigma_idx] * sigma_e[sigma_idx];
  }

  return V;
}
