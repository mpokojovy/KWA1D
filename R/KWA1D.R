#' Univariate Kernel Weighted Average Estimation
#'
#' Computes location and scale estimates for a univariate Student's t dataset via Kernel Weighted Average estimation
#'
#'
#'
#' @usage KWA1D.loc(x, df, bw = NULL, adjust = 1.0)
#' KWA1D.scale(x, df, bw.loc = NULL, adjust.loc = 1.0, bw.scale = NULL, adjust.scale = 1.0, bias.corr = NULL)
#' KWA1D.loc.scale(x, df, bw.loc = NULL, adjust.loc = 1.0, bw.scale = NULL, adjust.scale = 1.0, bias.corr = NULL)
#'
#' @param x numeric vector containing the dataset.
#' @param df degrees of freedom.
#' @param bw,bw.loc bandwidth for location estimation.
#' @param bw.scale bandwidth for scale estimation.
#' @param adjust,adjust.loc scalar adjustment to location bandwidth; default value is 1.0.
#' @param adjust.scale scalar adjustment to scale bandwidth; default value is 1.0.
#' @param bias.corr scalar factor to make the scale estimate unbiased.
#'
#' @details Consider a location-scale family \eqn{f(x|\mu, \sigma) = 1/\sigma \times f_0((x-\mu)/\sigma)} for real \eqn{x},
#' where \eqn{\mu} and \eqn{\sigma} are the location and scale parameters, respectively, and \eqn{f_0(\cdot)} is symmetric about 0.
#' For a squared-integrable random variable \eqn{X}, the location and scale parameters can be expressed as
#' \deqn{\mu = (\int x f^2(x) dx) /  (\int f^2(x) dx), }
#' \deqn{\sigma = ( ( \int z^2 f^2_0(z) dz) /  ( \int f^2_0 (z) dz) )^{-1/2} \times  ( ( \int (x-\mu)^2 f^2(x) dx) /  ( \int f^2(x) dx) )^{1/2}.}
#' These functionals can be estimated using appropriate KDE-based U-statistics leading to KWA estimates for location and scale parameters.
#'
#'
#' @return A scalar representing KWA location or scale estimate. When \code{KWA1D.loc.scale} is called, a list containing
#' location and scale estimates is returned.
#'
#' @author Michael Pokojovy (\email{michael.pokojovy@@gmail.com}), Su Chen, Andrews T. Anum, John Koomson
#'
#' @keywords Student's \emph{t}-distribution, kernel weighted average, KWA, bandwidth, kernel density functional, heavy tails
#'
#' @seealso \code{\link{KWA1D.estimator}}, \code{\link[robustbase]{Qn}}
#'
#' @references Pokojovy, M., Chen, S., Anum, A.T and Koomson, J. (2023) Adaptive location and scale estimation with kernel weighted averages. (submitted)
#'
#' Rousseeuw, P.J. and Croux, C. (1993) Alternatives to the median absolute deviation. \emph{Journal of the
#' American Statistical Association} \strong{88}, 1273–1283.
#'
#' Rousseeuw, P. J. and Leroy, A. M. (1987) \emph{Robust Regression and Outlier Detection}. Wiley.
#'
#'
#' @note If \code{bw} is not provided, \code{KWA1D.loc} finds an optimal bandwidth based on \code{df} and the sample size.
#' Similarly, \code{KWA1D.scale} and \code{KWA1D.loc.scale} compute \code{bw.loc} and \code{bw.scale} if they are not provided.
#'
#'
#'
# @source sources here
#'
#' @name KWA1D
#'
#'
#' @examples
#' n <- 100
#' df <- 2
#' x <- rt(n, df)
#' loc <- KWA1D.loc(x, df)
#' scale <- KWA1D.scale(x, df)
#' KWA1D.loc.scale(x, df)
#'
NULL

## Q_n scale estimator
.Qn_scale <- function(x, df) {
  constant = 1/(2.0^(0.5 + 0.5/df)*qt(5/8, df = df))
  return(robustbase::Qn(x, constant = constant, finite.corr = FALSE))
}

## Optimal bandwidth
KWA.opt.bw <- function(n, df, loc = TRUE) {
  df.grid = c(1:5, 10, 20, 30)

  if (df > max(df.grid)) {
    #warning(paste("df capped at", max(df.grid)))
    df = max(df.grid)
  }

  if (df < 1.0) {
    stop("df must be no less than 1")
  }

  if (loc) {
    coeff.opt = matrix(c(1.003823,  -0.1837867,
                         0.8743523,  0.000000,
                         1.245383,   0.000000,
                         1.52346,    0.000000,
                         1.748958,   0.000000,
                         2.387295,   0.000000,
                         3.155471,   0.000000,
                         3.591352,   0.000000), byrow = TRUE, ncol = 2)
  } else {
    coeff.opt = cbind(c(2.058333, 1.713362, 2.263782, 2.397197, 2.589637, 2.996948, 3.832673, 4.828874), 0.0)
  }

  ind = min(which(df.grid >= df))

  if (ind == 1) {
    h = coeff.opt[1, 1]*n^(coeff.opt[1, 2])
  } else {
    h1 = coeff.opt[ind - 1, 1]*n^(coeff.opt[ind - 1, 2])
    h2 = coeff.opt[ind,     1]*n^(coeff.opt[ind,     2])

    t = (df - df.grid[ind - 1])/(df.grid[ind] - df.grid[ind - 1])

    h = (1 - t)*h1 + t*h2
  }

  return(h)
}


## KWA scale bias correction
KWA.scale.bias.corr <- function(n, df) {
  df.grid = c(1:5, 10, 20, 30)

  if (df > max(df.grid)) {
    #warning(paste("df capped at", max(df.grid)))
    df = max(df.grid)
  }

  if (df < 1.0) {
    stop("df must be no less than 1")
  }

  coeff.opt = matrix(c(0.3906801,  -5.130337,   173.7469, -39.8597,  113.4416,
                       0.1207759,  -0.02644408, 1,         13.17681, 43.37158,
                       0.129282,    1.61891,    700.9076, -1.976418, 1,
                       0.1003642,  -0.3912979,  1,         21.71842, 65.27003,
                       0.08477952, -0.4686444,  1,         42.96931, 143.3675,
                       0.03365308, -0.8202625,  1,         30.75186, 86.17663,
                       0.01340489, -0.9055542,  1,         23.65226, 118.048,
                       0.01122529, -1.061898,   1,         23.51469, 77.98368), byrow = TRUE, ncol = 5)

  ind = min(which(df.grid >= df))

  if (ind == 1) {
    bias.corr = coeff.opt[1, 1] + coeff.opt[1, 2]/(n + coeff.opt[1, 3]) + coeff.opt[1, 4]/(n^2 + coeff.opt[1, 5]^2)
  } else {
    bias.corr1 = coeff.opt[ind - 1, 1] + coeff.opt[ind - 1, 2]/(n + coeff.opt[ind - 1, 3]) + coeff.opt[ind - 1, 4]/(n^2 + coeff.opt[ind - 1, 5]^2)
    bias.corr2 = coeff.opt[ind,     1] + coeff.opt[ind,     2]/(n + coeff.opt[ind,     3]) + coeff.opt[ind,     4]/(n^2 + coeff.opt[ind,     5]^2)

    t = (df - df.grid[ind - 1])/(df.grid[ind] - df.grid[ind - 1])

    bias.corr = (1 - t)*bias.corr1 + t*bias.corr2
  }

  bias.corr = exp(bias.corr)

  return(bias.corr)
}



## KWA loc
KWA1D.loc <- function(x, df, bw = NULL, adjust = 1.0) {
  if (is.null(df)) {
    stop("Provide df or use KWA1D.estimator() instead.")
  }

  if (df < 1.0) {
    stop("df must be no less than 1")
  }

  if (df <= 30) {
    if (is.null(bw)) {
      n = length(x)
      sigma.hat = .Qn_scale(x, df)
      bw = KWA.opt.bw(n = n, df = df, loc = TRUE)*sigma.hat
    }

    bw = adjust*bw

    return (KWA_loc_(x, bw))
  } else {
    return(mean(x))
  }
}

## KWA scale
KWA1D.scale <- function(x, df, bw.loc   = NULL, adjust.loc   = 1.0,
                        bw.scale = NULL, adjust.scale = 1.0,
                        bias.corr = KWA.scale.bias.corr(length(x), df)) {
  if (is.null(df)) {
    stop("Provide df or use KWA1D.estimator() instead.")
  }

  if (df < 1.0) {
    stop("df must be no less than 1")
  }

  if (length(x) < 2) {
    return(NA)
  }

  if (df <= 30) {
    if (is.null(bw.loc) || is.null(bw.scale)) {
      n = length(x)
      sigma.hat = .Qn_scale(x, df)
    }

    if (is.null(bw.loc)) {
      bw.loc = KWA.opt.bw(n = n, df = df, loc = TRUE)*sigma.hat
    }

    bw.loc = adjust.loc*bw.loc

    if (is.null(bw.scale)) {
      bw.scale = KWA.opt.bw(n = n, df = df, loc = FALSE)*sigma.hat
    }

    bw.scale = adjust.scale*bw.scale

    return(KWA_scale_(x, df, bw.loc, bw.scale)/bias.corr)
  } else {
    cf = if (df == Inf) 1.0 else sqrt((df - 2.0)/df)
    return(cf*sd(x))
  }
}

## KWA loc and scale
KWA1D.loc.scale <- function(x, df, bw.loc = NULL, adjust.loc   = 1.0,
                            bw.scale = NULL, adjust.scale = 1.0,
                            bias.corr = KWA.scale.bias.corr(length(x), df)) {
  if (is.null(df)) {
    stop("Provide df or use KWA1D.estimator() instead.")
  }

  if (df <= 30) {
    n = length(x)

    if (is.null(bw.loc) || is.null(bw.scale)) {
      n = length(x)
      sigma.hat = .Qn_scale(x, df)
    }

    if (is.null(bw.loc)) {
      bw.loc = KWA.opt.bw(n = n, df = df, loc = TRUE)*sigma.hat
    }

    bw.loc = adjust.loc*bw.loc

    if (is.null(bw.scale)) {
      bw.scale = KWA.opt.bw(n = n, df = df, loc = FALSE)*sigma.hat
    }

    bw.scale = adjust.scale*bw.scale

    return(list(loc   = KWA_loc_(x, bw.loc),
                scale = KWA_scale_(x, df, bw.loc, bw.scale)/bias.corr))
  } else {
    cf = if (df == Inf) 1.0 else sqrt((df - 2.0)/df)
    return(list(loc = mean(x), scale = cf*sd(x)))
  }
}
