

#'  Model R-squared
#'
#' computes the R-squared for several generalized linear models

#' @param object object of class "glm", or an object for which the function \code{\link[stats]{deviance}}
#'          or \code{\link[stats]{deviance}} is defined.
#' @param test logical. If TRUE, include the associated ANOVA table.
#' @param ci logical. If TRUE, compute a confidence interval for the primary index.
#' @param ci_width numeric confidence level between 0 and 1.
#' @param quiet logical. If TRUE, muffle warnings and messages emitted while
#'   computing the index.
#' @param ... not implemented yet
#' @return A list with `indices`, optional `anova` and `ci` tables, and the
#'   deviance or sum-of-squares quantities used in the calculations.
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' r2(model)
#'
#' @rdname r2
#' @export

r2 <- function(object, test = FALSE, ci = FALSE, ci_width = 0.95,
               quiet = FALSE, ...) {
  if (isTRUE(quiet)) {
    return(withCallingHandlers(
      r2(object, test=test, ci=ci, ci_width=ci_width, quiet=FALSE, ...),
      warning=function(e) invokeRestart("muffleWarning"),
      message=function(e) invokeRestart("muffleMessage")
    ))
  }
  UseMethod("r2")
}

#' @rdname r2
#' @export

r2.default<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                     quiet=FALSE,...) {


  d1<-stats::deviance(object)
  if (is.null(d1)) {
    d1<- -2*stats::logLik(object)
  }
  model0<-stats::update(object,.~1,data=object$model)
  df<-length(stats::coef(object))-length(stats::coef(model0))
  d0<-stats::deviance(model0)
  if (is.null(d0))
    d0<- -2*stats::logLik(model0)

  r2<-1-d1/d0
  r2adj<-1-(d1+df)/d0
  obj<-data.frame(r2=r2,r2adj=r2adj)
  rownames(obj)<-"Model"

  tab<-NULL
  if (isTRUE(test)) {
    # Compare the reduced model to the fitted model so the LR difference is
    # reported as null deviance minus residual deviance.
    tab<-stats::anova(model0,object)
    if (nrow(tab)>1)  tab<-tab[2,,drop=FALSE]
    names(tab)<-transnames(names(tab),list(test=c("Deviance","LR stat.","LR.stat"),p=c("Pr(>Chi)","Pr(Chi)","Pr(>Chisq)")))
  }
  ci_tab <- NULL
  if (isTRUE(ci)) {
    ci_tab <- .ci_eta2_table(obj$r2, df, d0, ci_width,
                             row_names=rownames(obj))
    names(ci_tab)[1] <- "r2"
  }
  structure(list(
    indices=obj,
    anova=tab,
    ci=ci_tab,
    df=df,
    D0=d0,
    D1=d1,
    D=d1,
    model="glm"
  ),class="nwa"
  )

}

#' @rdname r2
#' @export

r2.lm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                quiet=FALSE,...) {

  ss<-summary(object)
  df1 <- unname(ss$fstatistic["numdf"])
  mf <- stats::model.frame(object)
  y <- stats::model.response(mf)
  weights <- stats::model.weights(mf)
  if (is.null(weights))
    weights <- rep(1, length(y))
  ybar <- stats::weighted.mean(y, weights)
  d0 <- sum(weights * (y - ybar)^2)
  d1 <- stats::deviance(object)
  obj<-data.frame(r2=  ss$r.squared,r2adj=  ss$adj.r.squared)
  rownames(obj)<-"Model"
  tab<-NULL
  if (isTRUE(test)) {
    f <- unname(ss$fstatistic["value"])
    df2 <- unname(ss$fstatistic["dendf"])
    p <- stats::pf(
      f,
      df1 = df1,
      df2 = df2,
      lower.tail = FALSE
    )
    tab<-data.frame(f=f,df=df1,dfr=df2,p=p)
  }
  ci_tab <- NULL
  if (isTRUE(ci)) {
    ci_tab <- .ci_eta2_table(obj$r2, df1, d0, ci_width,
                             row_names=rownames(obj))
    names(ci_tab)[1] <- "r2"
  }
  structure(list(
    indices=obj,
    anova=tab,
    ci=ci_tab,
    df=df1,
    D0=d0,
    D1=d1,
    D=d1,
    model="lm"
  ),
  class="nwa"
  )
}

#' @rdname r2
#' @export

r2.glm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                 quiet=FALSE,...) {

  r2.default(object,test=test,ci=ci,ci_width=ci_width,quiet=quiet,...)
}

#' @rdname r2
#' @export


r2.clm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                 quiet=FALSE,...) {

  r2.default(object,test=test,ci=ci,ci_width=ci_width,quiet=quiet,...)
}

#' @rdname r2
#' @export

r2.multinom<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                      quiet=FALSE,...) {

  if (is.null(object$model))
    stop("model of class `multinom` should be estimated with `nnet::multinom(...,model=TRUE)` option")
  r2.default(object,test=test,ci=ci,ci_width=ci_width,quiet=quiet,...)
}

#'  Eta-squared and Epsilon-squared
#'
#' Computes the eta-squared and epsilon-squared indices for several generalized linear models

#' @param object object of class "glm", or an object for which the function \code{\link[stats]{drop1}}
#'          or \code{\link[stats]{deviance}} is defined.
#' @param test logical. If TRUE, include the associated ANOVA table.
#' @param ci logical. If TRUE, compute confidence intervals for the indices.
#' @param ci_width numeric confidence level between 0 and 1.
#' @param quiet logical. If TRUE, muffle warnings and messages emitted while
#'   computing the indices.
#' @param .test character. Test statistic passed to `car::Anova()` where applicable.
#' @param col character. Column of the ANOVA table used as the effect statistic.
#' @param anova_table optional precomputed ANOVA table.
#' @param ... not implemented yet
#' @return A list with `indices`, optional `anova` and `ci` tables, and the
#'   deviance or sum-of-squares quantities used in the calculations.
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' eta2(model)
#' @rdname eta2
#' @export

eta2 <- function(object, test = FALSE, ci = FALSE, ci_width = 0.95,
                 quiet = FALSE, ...) {
  if (isTRUE(quiet)) {
    return(withCallingHandlers(
      eta2(object, test=test, ci=ci, ci_width=ci_width, quiet=FALSE, ...),
      warning=function(e) invokeRestart("muffleWarning"),
      message=function(e) invokeRestart("muffleMessage")
    ))
  }
  UseMethod("eta2")
}

#' @rdname eta2
#' @export

eta2.default<-function(object, test=FALSE, ci=FALSE, ci_width=0.95,
                       quiet=FALSE, .test="LR", col=NULL,
                       anova_table=NULL, ...) {

  if (is.null(col)) {
    col <- switch(.test,
      Wald = "Chisq",
      F = "F",
      Chisq = "Chisq",
      "LR Chisq"
    )
  }

  model0<-stats::update(object ,.~1,data=object$model)
  dev0<-stats::deviance(model0)
  if (is.null(dev0))
     dev0<- as.numeric(-2*stats::logLik(model0))
  a <- anova_table
  if (is.null(a)) {
    .test <- match.arg(.test, c("LR", "Wald", "F"))
    a<-car::Anova(object,type=3,test=.test)
  }
  df<-a$Df
  # etas
  res<-matrix(a[,col]/dev0,ncol = 1)
  rownames(res)<-rownames(a)
  colnames(res)<-"Eta_squared"
  #gammas
  gam<-matrix((a[,col]-df)/dev0,ncol = 1)
  rownames(gam)<-rownames(a)
  colnames(gam)<-"Epsilon_squared"
  gam[gam < 0] <- 0
  obj<-cbind(res,gam)
  ci_tab <- NULL
  if (isTRUE(ci)) {
    ci_tab <- .ci_eta2_table(obj[, "Eta_squared"], df, dev0, ci_width,
                             row_names=rownames(obj))
  }

  structure(list(
    indices=obj,
    anova=if (isTRUE(test)) a else NULL,
    ci=ci_tab,
    df=df,
    D0=dev0,
    D1=a[,col],
    D=a[,col],
    model="glm"
  ),
  class="nwa"
  )


}

#' @rdname eta2
#' @export

eta2.lm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                  quiet=FALSE, anova_table=NULL,...) {

  a<-if (is.null(anova_table)) car::Anova(object,type=3) else anova_table
  model0<-stats::update(object ,.~1,data=object$model)
  w<-which(rownames(a) %in% c("(Intercept)","Residuals"))

  sse0<-model0$df.residual*stats::sigma(model0)^2
  msem<-a$`Sum Sq`[nrow(a)]/object$df.residual   # Residuals row is last; MSE of the full model

  res<-a$`Sum Sq`/sse0
  res<-matrix(res[-w],ncol=1)
  rownames(res)<-rownames(a)[-w]
  colnames(res)<-"Eta_squared"

  eps<-(a$`Sum Sq`-a$Df*msem)/sse0
  eps<-matrix(eps[-w],ncol=1)
  rownames(eps)<-rownames(a)[-w]
  colnames(eps)<-"Epsilon_squared"
  eps[eps < 0] <- 0
  effect_ss <- a$`Sum Sq`[-w]
  obj<-cbind(res,eps)
  ci_tab <- NULL
  if (isTRUE(ci)) {
    ci_tab <- .ci_eta2_table(obj[, "Eta_squared"], a$Df[-w], sse0,
                             ci_width, row_names=rownames(obj))
  }

  structure(list(
    indices=obj,
    anova=if (isTRUE(test)) a else NULL,
    ci=ci_tab,
    df=a$Df[-w],
    D0=sse0,
    D1=effect_ss,
    D=effect_ss,
    model="lm"
  ),
  class="nwa"
  )


}

#' @rdname eta2
#' @export

eta2.glm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                   quiet=FALSE, .test="LR",col=NULL,
                   anova_table=NULL,...) {

  eta2.default(object,test=test,ci=ci,ci_width=ci_width,
               quiet=quiet,.test=.test,col=col,
               anova_table=anova_table,...)

}

#' @rdname eta2
#' @export

eta2.clm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                   quiet=FALSE, col="Chisq",anova_table=NULL,...) {

  # car::Anova()/ordinal::clm's own anova() both compute a Wald test for a single `clm`
  # fit (there is no refit-based Type III LR test available for this class), so a genuine
  # per-term likelihood-ratio table is built directly instead; see .clm_anova_lr().
  if (is.null(anova_table))
    anova_table <- .clm_anova_lr(object)
  eta2.default(object,test=test,ci=ci,ci_width=ci_width,
               quiet=quiet,col=col,anova_table=anova_table,...)

}

#' @rdname eta2
#' @export

eta2.multinom<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                        quiet=FALSE, .test="LR",col="LR Chisq",
                        anova_table=NULL,...) {

  if (is.null(object$model))
    stop("model of class `multinom` should be estimated with `nnet::multinom(...,model=TRUE)` option")
  eta2.default(object,test=test,ci=ci,ci_width=ci_width,
               quiet=quiet,.test=.test,col=col,
               anova_table=anova_table,...)

}



#'  Partial Eta-squared and Partial Epsilon-squared
#'
#' Computes the partial eta-squared and partial epsilon-squared indices for several generalized linear models

#' @param object object of class "glm", or an object for which the function \code{\link[stats]{drop1}}
#'          or \code{\link[stats]{deviance}} is defined.
#' @param test logical. If TRUE, include the associated ANOVA table.
#' @param ci logical. If TRUE, compute confidence intervals for the indices.
#' @param ci_width numeric confidence level between 0 and 1.
#' @param quiet logical. If TRUE, muffle warnings and messages emitted while
#'   computing the indices.
#' @param .test character. Test statistic passed to `car::Anova()` where applicable.
#' @param col character. Column of the ANOVA table used as the effect statistic.
#' @param anova_table optional precomputed ANOVA table.
#' @param ... not implemented yet
#' @return A list with `indices`, optional `anova` and `ci` tables, and the
#'   deviance or sum-of-squares quantities used in the calculations.
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' eta2_partial(model)
#' @rdname eta2_partial
#' @export

eta2_partial <- function(object, test = FALSE, ci = FALSE,
                         ci_width = 0.95, quiet = FALSE, ...) {
  if (isTRUE(quiet)) {
    return(withCallingHandlers(
      eta2_partial(object, test=test, ci=ci, ci_width=ci_width,
                   quiet=FALSE, ...),
      warning=function(e) invokeRestart("muffleWarning"),
      message=function(e) invokeRestart("muffleMessage")
    ))
  }
  UseMethod("eta2_partial")
}

#' @rdname eta2_partial
#' @export

eta2_partial.default<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                                quiet=FALSE, .test="LR",col=NULL,
                                anova_table=NULL,...) {

  if (is.null(col)) {
    col <- switch(.test,
      Wald = "Chisq",
      F = "F",
      Chisq = "Chisq",
      "LR Chisq"
    )
  }

  devm<-stats::deviance(object)
  if (is.null(devm))
    devm<- as.numeric(-2*stats::logLik(object))
  a <- anova_table
  if (is.null(a)) {
    .test <- match.arg(.test, c("LR", "Wald", "F"))
    a<-car::Anova(object,type=3,test=.test)
  }
  df<-a$Df
  k<-sum(df)
  # D_{m.x}, the deviance of the model without each term, obtained from the
  # term's LR chi-squared (D_{m.x} - D_m) plus the full model deviance D_m
  devmx<-devm+a[,col]
  # petas: (D_{m.x}-D_m)/D_{m.x}
  res<-matrix(a[,col]/devmx,ncol = 1)
  rownames(res)<-rownames(a)
  colnames(res)<-"Eta2_p"
  #gammas: (D_{m.x}-D_m-u)/(D_{m.x}+k-u)
  gam<-matrix((a[,col]-df)/(devmx+k-df),ncol = 1)
  rownames(gam)<-rownames(a)
  colnames(gam)<-"Epsilon2_p"
  gam[gam < 0] <- 0
  obj<-cbind(res,gam)
  ci_tab <- NULL
  if (isTRUE(ci)) {
    ci_tab <- .ci_eta2p_table(obj[, "Eta2_p"], df, devmx, ci_width,
                              row_names=rownames(obj))
  }
  structure(list(
    indices=obj,
    anova=if (isTRUE(test)) a else NULL,
    ci=ci_tab,
    df=df,
    D0=devm,
    D1=a[[col]],
    D=a[[col]],
    model="glm"
  ),
  class="nwa"
  )
}

#' @rdname eta2_partial
#' @export

eta2_partial.lm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                           quiet=FALSE, anova_table=NULL,...) {

  a <- if (is.null(anova_table)) car::Anova(object,type=3) else anova_table
  w <- which(rownames(a) %in% c("(Intercept)","Residuals"))
  effect_ss <- a$`Sum Sq`[-w]
  effect_df <- a$Df[-w]
  sse_full <- stats::deviance(object)
  devmx <- sse_full + effect_ss
  k <- sum(effect_df)
  msem <- sse_full / object$df.residual

  res <- matrix(effect_ss / devmx, ncol=1)
  rownames(res) <- rownames(a)[-w]
  colnames(res) <- "Eta2_p"

  eps <- matrix(
    (effect_ss - effect_df * msem) / (devmx + k - effect_df),
    ncol=1
  )
  rownames(eps) <- rownames(a)[-w]
  colnames(eps) <- "Epsilon2_p"
  eps[eps < 0] <- 0
  obj <- cbind(res, eps)
  ci_tab <- NULL
  if (isTRUE(ci)) {
    ci_tab <- .ci_eta2p_table(obj[, "Eta2_p"], effect_df, devmx,
                              ci_width, row_names=rownames(obj))
  }

  structure(list(
    indices=obj,
    anova=if (isTRUE(test)) a else NULL,
    ci=ci_tab,
    df=effect_df,
    D0=sse_full,
    D1=effect_ss,
    D=effect_ss,
    model="lm"
  ),
  class="nwa"
  )
}

#' @rdname eta2_partial
#' @export

eta2_partial.glm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                            quiet=FALSE, .test="LR",col=NULL,
                            anova_table=NULL,...) {

  eta2_partial.default(object,test=test,ci=ci,ci_width=ci_width,
                       quiet=quiet,.test=.test,col=col,
                       anova_table=anova_table,...)
}

#' @rdname eta2_partial
#' @export

eta2_partial.clm<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                            quiet=FALSE, col="Chisq",
                            anova_table=NULL,...) {

  # see eta2.clm(): car::Anova()/clm's own anova() are Wald-based for a single fit, so a
  # genuine per-term likelihood-ratio table is built directly via .clm_anova_lr().
  if (is.null(anova_table))
    anova_table <- .clm_anova_lr(object)
  eta2_partial.default(object,test=test,ci=ci,ci_width=ci_width,
                        quiet=quiet,col=col,anova_table=anova_table,...)

}

#' @rdname eta2_partial
#' @export

eta2_partial.multinom<-function(object,test=FALSE,ci=FALSE,ci_width=0.95,
                                 quiet=FALSE, .test="LR",col="LR Chisq",
                                 anova_table=NULL,...) {

  if (is.null(object$model))
    stop("model of class `multinom` should be estimated with `nnet::multinom(...,model=TRUE)` option")
  eta2_partial.default(object,test=test,ci=ci,ci_width=ci_width,
                       quiet=quiet,.test=.test,col=col,
                       anova_table=anova_table,...)

}

#'  Print numeric with attributes
#'
#'  Prints generic named vectors without printing their attributes

#' @param x object of class "nwa" (numeric with attributes)
#'          or \code{\link[stats]{deviance}} is defined.
#' @param ... not implemented yet
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' print(r2(model))
#' @export

print.nwa<-function(x,...) {
 indices <- x$indices
 if (!is.null(x$ci)) {
   ci <- x$ci
   if (!is.null(rownames(indices)) && !is.null(rownames(ci)))
     ci <- ci[rownames(indices), , drop=FALSE]
   ci_side <- ci[, c("lower", "upper"), drop=FALSE]
   names(ci_side) <- c("CI_lower", "CI_upper")
   indices <- cbind(
     indices[, 1, drop=FALSE],
     ci_side,
     indices[, -1, drop=FALSE]
   )
 }
 print(indices)
 if (!is.null(x$anova)) {
   cat("\n")
   print(x$anova)
 }
}
