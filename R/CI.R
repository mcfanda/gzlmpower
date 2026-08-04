#' Confidence interval for eta-squared / via noncentral chi-square inversion
#'
#' @param eta2  numeric. the effect size index for predictor x
#' @param u     integer. Degrees of freedom associated with predictor x.
#' @param D0    numeric. Null-model deviance. Required for eta2;
#' @param conf.level numeric. Confidence level (default 0.95).
#'
#' @return A data.frame  with eta2, lower and upper confidence limits.
#'
#' @examples
#' ci_eta2(eta2 = .20, u = 2, D0 = 1000)
#' @export

ci_eta2 <- function(eta2, u, D0 , conf.level = 0.95) {

  if (!requireNamespace("MBESS", quietly = TRUE)) {
    stop("Package 'MBESS' is required. Install it with install.packages('MBESS').")
  }
  if (is.null(D0)) {
    stop("Supply  D0 (for eta2).")
  }
  Qx<-eta2*D0
  # Noncentrality-parameter CI via inversion of the noncentral chi-square CDF.
  # MBESS::conf.limits.nc.chisq() clips lambda.L at 0 internally when Qx < u,
  # matching the boundary rule described in the manuscript.
  nc_ci <- MBESS::conf.limits.nc.chisq(
    Chi.Square = Qx,
    df         = u,
    conf.level = conf.level
  )

  lambda_L <- nc_ci$Lower.Limit
  lambda_U <- nc_ci$Upper.Limit

  # NA handling: conf.limits.nc.chisq() can return NULL/NA at the boundary;
  # treat that as lambda_L = 0 per the manuscript's convention.
  if (is.null(lambda_L) || is.na(lambda_L)) lambda_L <- 0

  results<-data.frame(eta2=eta2,lower=lambda_L / D0,upper = lambda_U / D0)

  attr(results,"Qx")<-Qx
  attr(results,"u")<-u
  attr(results,"lambda.L")<-lambda_L
  attr(results,"lambda.U")<-lambda_U
  attr(results,"conf.level")<-conf.level
  results
}


#' Confidence interval for partial eta-squared
#'
#' @param eta2p  numeric. the effect size index for predictor x
#' @param u     integer. Degrees of freedom associated with predictor x.
#' @param Dmx   numeric. Deviance of the model excluding predictor x.
#' @param conf.level numeric. Confidence level (default 0.95).
#'
#' @return A data.frame  with eta2, lower and upper confidence limits.
#'
#' @examples
#' ci_eta2p(eta2p = .20, u = 2, Dmx = 1000)
#' @export

ci_eta2p <- function(eta2p, u, Dmx , conf.level = 0.95) {

  if (!requireNamespace("MBESS", quietly = TRUE)) {
    stop("Package 'MBESS' is required. Install it with install.packages('MBESS').")
  }
  if (is.null(Dmx)) {
    stop("Supply Dmx for partial eta2.")
  }
  Qx<-eta2p*Dmx
  # Noncentrality-parameter CI via inversion of the noncentral chi-square CDF.
  # MBESS::conf.limits.nc.chisq() clips lambda.L at 0 internally when Qx < u,
  # matching the boundary rule described in the manuscript.
  nc_ci <- MBESS::conf.limits.nc.chisq(
    Chi.Square = Qx,
    df         = u,
    conf.level = conf.level
  )

  lambda_L <- nc_ci$Lower.Limit
  lambda_U <- nc_ci$Upper.Limit

  # NA handling: conf.limits.nc.chisq() can return NULL/NA at the boundary;
  # treat that as lambda_L = 0 per the manuscript's convention.
  if (is.null(lambda_L) || is.na(lambda_L)) lambda_L <- 0

  results<-data.frame(eta2p=eta2p,lower=lambda_L / Dmx,upper = lambda_U / Dmx)

  attr(results,"Qx")<-Qx
  attr(results,"u")<-u
  attr(results,"lambda.L")<-lambda_L
  attr(results,"lambda.U")<-lambda_U
  attr(results,"conf.level")<-conf.level
  results
}



# CI helpers are evaluated only when the user requests confidence intervals.
.validate_ci_width <- function(ci_width) {
  if (!is.numeric(ci_width) || length(ci_width) != 1L ||
      !is.finite(ci_width) || ci_width <= 0 || ci_width >= 1) {
    stop("ci_width must be a single number between 0 and 1")
  }
  ci_width
}

.ci_eta2_table <- function(estimate, df, D0, conf.level, row_names = NULL) {
  conf.level <- .validate_ci_width(conf.level)
  estimate <- as.numeric(estimate)
  df <- rep(df, length.out = length(estimate))
  D0 <- rep(D0, length.out = length(estimate))
  out <- do.call(rbind, lapply(seq_along(estimate), function(i) {
    ci_eta2(estimate[i], df[i], D0[i], conf.level)
  }))
  if (is.null(row_names))
    row_names <- as.character(seq_along(estimate))
  rownames(out) <- row_names
  out
}

.ci_eta2p_table <- function(estimate, df, Dmx, conf.level, row_names = NULL) {
  conf.level <- .validate_ci_width(conf.level)
  estimate <- as.numeric(estimate)
  df <- rep(df, length.out = length(estimate))
  Dmx <- rep(Dmx, length.out = length(estimate))
  out <- do.call(rbind, lapply(seq_along(estimate), function(i) {
    ci_eta2p(estimate[i], df[i], Dmx[i], conf.level)
  }))
  if (is.null(row_names))
    row_names <- as.character(seq_along(estimate))
  rownames(out) <- row_names
  out
}
