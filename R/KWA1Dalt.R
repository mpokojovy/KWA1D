#' Adaptive Univariate Kernel Weighted Average Estimation
#'
#' Computes Kernel Weighted Average estimates of univariate Student's \emph{t} location and scale
#'
#' @usage KWA1D.estimator(x, df.search)
#'
# @inheritParams KWA1D
#' @param x numeric vector containing the dataset.
#' @param df.search numeric vector representing a grid to search for ``best" \code{df} match.
#'
#' @details See \code{\link{KWA1D}}
#'
#' @return  A list consisting of KWA location and scale estimates together with the ``best" \code{df} and associated maximal
#' \code{p-value} from the Cramer-von Mises goodness-of-fit test.
#'
#' @author Michael Pokojovy (\email{michael.pokojovy@@gmail.com}), Su Chen, Andrews T. Anum, John Koomson
#' @keywords Student's \emph{t}-distribution, kernel weighted average, KWA, bandwidth, kernel density functional, heavy tails
#'
#' @seealso \code{\link{KWA1D}}, \code{\link[robustbase]{Qn}}
#'
#' @references Pokojovy, M., Chen, S., Anum, A.T and Koomson, J. (2023) Adaptive location and scale estimation with kernel weighted averages. (submitted)
#'
#' Rousseeuw, P.J. and Croux, C. (1993) Alternatives to the Median Absolute Deviation. \emph{Journal of the
#' American Statistical Association} \strong{88}, 1273–1283.
#'
#' Rousseeuw, P. J. and Leroy, A. M. (1987) \emph{Robust Regression and Outlier Detection}. Wiley.
#'
#' @note If \code{df.search} is not provided, the default search interval is
#' \code{\strong{c(seq(from = 1.0, to = 5,  by = 0.2), seq(from = 5.5, to = 15, by = 0.5), seq(from = 16,  to = 30, by = 1),
#' 50, 100, 500, Inf)}}
#'
#'
#' @name KWA1D.estimator
#'
#' @examples
#' n <- 100
#' df <- 2
#' x <- rt(n, df)
#' KWA1D.estimator(x) #using default search interval
#' KWA1D.estimator(x, df.search = seq(1, 5, by = 0.05)) #using a user-provided search interval
#'

## Estimate df using Cramer-von Mises test
KWA1D.estimator <- function(x, df.search = c(seq(from = 1.0, to = 5,  by = 0.2),
                                             seq(from = 5.5, to = 15, by = 0.5),
                                             seq(from = 16,  to = 30, by = 1), 50, 100, 500, Inf)) {
  loc   = NA
  scale = NA
  pval  = 0.0
  df    = NA

  if (min(df.search) < 1) {
    stop("df.search cannot contain values less than 1.0")
  }

  for (i.df in 1:length(df.search)) {
    .df = df.search[i.df]

    est = KWA1D.loc.scale(x, .df)

    z = (x - est$loc)/est$scale

    .pval = goftest::cvm.test(z, "pt", df = .df)$p.value

    if (.pval >= pval) {
      pval  = .pval
      df    = .df
      loc   = est$loc
      scale = est$scale
    }
  }

  return(list(loc = loc, scale = scale, df = df, pval = pval))
}
