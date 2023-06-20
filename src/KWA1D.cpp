#include <Rcpp.h>
using namespace Rcpp;

static double C = 1.0/::sqrt(2*M_PI);

double phi(double x) {
  return C*::exp(-0.5*x*x);
}

// [[Rcpp::export]]
double KWA_loc_(const Rcpp::NumericVector &x, double bw) {
  double num = 0.0;
  double den = 0.0;

  const unsigned int n = x.size();
  const double c = 0.5/(bw*n*(n - 1.0));

  for (unsigned int i = 0; i < n; i++) {
    const double xi = x[i];

    for (unsigned int j = 0; j < i; j++) {
      const double xj = x[j];
      const double K = phi((xi - xj)/bw);

      num += c*(xi + xj)*K;
      den += 2.0*c*K;
    }
  }

  return num/den;
}


// [[Rcpp::export]]
double KWA_scale_(const Rcpp::NumericVector x, double df, double bw_loc, double bw_scale) {
  const unsigned int n = x.size();

  const double loc = KWA_loc_(x, bw_loc);

  const double c = 1.0/(bw_scale*n*(n - 1.0));

  double num = 0.0;
  double den = 0.0;

  for (unsigned int i = 0; i < n; i++) {
    const double xi  = x[i];
    const double di = xi - loc;

    for (unsigned int j = 0; j < i; j++) {
      const double xj = x[j];

      const double K = phi((xi - xj)/bw_scale);

      const double dj = xj - loc;

      num += c*(di*di + dj*dj)*K;
      den += 2.0*c*K;
    }
  }

  return ::sqrt(std::max(0.0, num/den));
}
