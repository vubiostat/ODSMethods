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
  const int nrowzi = zi.nrow();
  const int ncolzi = zi.ncol();
  const int nERRsd = sigma_e.size();

  // Construct random-effects correlation matrix.
  //
  // R:
  // b <- matrix(0, ncolzi, ncolzi)
  // b[lower.tri(b)] <- rho.vc
  // CorMat.RE <- t(b) + b + diag(ncolzi)

  Rcpp::NumericMatrix CorMat_RE(ncolzi, ncolzi);

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

  // V = Z D Z'
  //
  // D[i,j] = sigma_vc[i] *
  //          CorMat_RE[i,j] *
  //          sigma_vc[j]
  //
  // Rather than constructing D, incorporate the standard deviations
  // directly into the matrix multiplication.

  Rcpp::NumericMatrix V(nrowzi, nrowzi);

  for (int i = 0; i < nrowzi; ++i)
  {
    for (int j = 0; j < nrowzi; ++j)
    {
      double value = 0.0;

      for (int k = 0; k < ncolzi; ++k)
      {
        for (int l = 0; l < ncolzi; ++l)
        {
          value +=
            zi(i, k) *
            sigma_vc[k] *
            CorMat_RE(k, l) *
            sigma_vc[l] *
            zi(j, l);
        }
      }

      V(i, j) = value;
    }
  }

  // Add measurement-error variance to the diagonal.
  //
  // R:
  // diag(rep(sigma.e^2, each = nrowzi / nERRsd))

  const int observations_per_sigma = nrowzi / nERRsd;

  for (int i = 0; i < nrowzi; ++i)
  {
    const int sigma_idx = i / observations_per_sigma;

    V(i, i) += sigma_e[sigma_idx] * sigma_e[sigma_idx];
  }

  return V;
}
