
#'  Model R-squared
#'
#' computes the R-squared for several generalized linear models

#' @param object object of class "glm", or an object for which the function \code{\link[stats]{deviance}}
#'          or \code{\link[stats]{deviance}} is defined.
#' @param ... not implemented yet
#' @return an numeric vector with the R-square and adjusted R-square
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' r2(model)
#'
#' @rdname r2
#' @export

r2 <- function(object, ...) UseMethod("r2")

#' @rdname r2
#' @export

r2.default<-function(object,test=FALSE,...) {


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
  if (test) {
    tab<-anova(object,model0)
    if (nrow(tab)>1)  tab<-tab[2,,drop=FALSE]
    names(tab)<-transnames(names(tab),list(test=c("Deviance","LR stat.","LR.stat"),p=c("Pr(>Chi)","Pr(Chi)","Pr(>Chisq)")))
    obj$test<-tab$test
    obj$p<-tab$p
  }
  attr(obj,"df")<-df
  class(obj)<-c(class(obj),"nwa")
  obj
}

#' @rdname r2
#' @export

r2.lm<-function(object,test=FALSE,...) {

  ss<-summary(object)
  res<-data.frame(r2=  ss$r.squared,r2adj=  ss$adj.r.squared)
  if (test) {
    sm<-summary(object)
    f <- unname(sm$fstatistic["value"])
    df1 <- unname(sm$fstatistic["numdf"])
    df2 <- unname(sm$fstatistic["dendf"])
    p <- stats::pf(
      f,
      df1 = df1,
      df2 = df2,
      lower.tail = FALSE
    )
    res$test<-f
    res$p<-p
  }
  res
}

#' @rdname r2
#' @export

r2.glm<-function(object,test=FALSE,...) {

  r2.default(object,test=test)
}

#' @rdname r2
#' @export


r2.clm<-function(object,test=FALSE,...) {

  r2.default(object,test=test)
}

#' @rdname r2
#' @export

r2.multinom<-function(object,test=FALSE,...) {

  if (is.null(object$model))
    stop("model of class `multinom` should be estimated with `nnet::multinom(...,model=TRUE)` option")
  r2.default(object,test=test)
}

#'  Eta-squared and Epsilon-squared
#'
#' Computes the eta-squared and epsilon-squared indices for several generalized linear models

#' @param object object of class "glm", or an object for which the function \code{\link[stats]{drop1}}
#'          or \code{\link[stats]{deviance}} is defined.
#' @param ... not implemented yet
#' @return an anova table with R
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' eta2(model)
#' @rdname eta2
#' @export

eta2 <- function(object, ...) UseMethod("eta2")

#' @rdname eta2
#' @export

eta2.default<-function(object,...) {

  args<-list(...)
  test<-"LR"
  col<-"LR Chisq"

  if (utils::hasName(args,"test"))
     test<-args$test
  if (utils::hasName(args,"col"))
    col<-args$col

  model0<-stats::update(object ,.~1,data=object$model)
  dev0<-stats::deviance(model0)
  if (is.null(dev0))
     dev0<- as.numeric(-2*stats::logLik(model0))
  a<-car::Anova(object,type=3,test=test)
  df<-a$Df
  # etas
  res<-matrix(a[,col]/dev0,ncol = 1)
  rownames(res)<-rownames(a)
  colnames(res)<-"Eta_squared"
  #gammas
  gam<-matrix((a[,col]-df)/dev0,ncol = 1)
  rownames(gam)<-rownames(a)
  colnames(gam)<-"Epsilon_squared"
  gam[gam<0]<-0
  cbind(res,gam)
}

#' @rdname eta2
#' @export

eta2.lm<-function(object,...) {

  a<-car::Anova(object,type=3)
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
  eps[eps<0]<-0
  cbind(res,eps)

}

#' @rdname eta2
#' @export

eta2.glm<-function(object,...) {

  eta2.default(object,test="LR")

}

#' @rdname eta2
#' @export

eta2.clm<-function(object,...) {

  eta2.default(object,test="Chisq",col="Chisq")

}

#' @rdname eta2
#' @export

eta2.multinom<-function(object,...) {

  if (is.null(object$model))
    stop("model of class `multinom` should be estimated with `nnet::multinom(...,model=TRUE)` option")
  eta2.default(object,test="Chisq",col="LR Chisq")

}



#'  Partial Eta-squared and Partial Epsilon-squared
#'
#' Computes the partial eta-squared and partial epsilon-squared indices for several generalized linear models

#' @param object object of class "glm", or an object for which the function \code{\link[stats]{drop1}}
#'          or \code{\link[stats]{deviance}} is defined.
#' @param ... not implemented yet
#' @return an anova table with R
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' eta2_partial(model)
#' @rdname eta2_partial
#' @export

eta2_partial <- function(object, ...) UseMethod("eta2_partial")

#' @rdname eta2
#' @export

eta2_partial.default<-function(object,...) {

  args<-list(...)
  test<-"LR"
  col<-"LR Chisq"

  if (utils::hasName(args,"test"))
    test<-args$test
  if (utils::hasName(args,"col"))
    col<-args$col

  devm<-stats::deviance(object)
  if (is.null(devm))
    devm<- as.numeric(-2*stats::logLik(object))
  a<-car::Anova(object,type=3,test=test)
  df<-a$Df
  k<-sum(df)
  # D_{m.x}, the deviance of the model without each term, obtained from the
  # term's LR chi-squared (D_{m.x} - D_m) plus the full model deviance D_m
  devmx<-devm+a[,col]
  # petas: (D_{m.x}-D_m)/D_{m.x}
  res<-matrix(a[,col]/devmx,ncol = 1)
  rownames(res)<-rownames(a)
  colnames(res)<-"Eta_squared"
  #gammas: (D_{m.x}-D_m-u)/(D_{m.x}+k-u)
  gam<-matrix((a[,col]-df)/(devmx+k-df),ncol = 1)
  rownames(gam)<-rownames(a)
  colnames(gam)<-"Epsilon_squared"
  gam[gam<0]<-0
  cbind(res,gam)
}

#' @rdname eta2_partial
#' @export

eta2_partial.clm<-function(object,...) {

  eta2_partial.default(object,test="Chisq",col="Chisq")

}

#' @rdname eta2_partial
#' @export

eta2_partial.multinom<-function(object,...) {

  if (is.null(object$model))
    stop("model of class `multinom` should be estimated with `nnet::multinom(...,model=TRUE)` option")
  eta2_partial.default(object,test="Chisq",col="LR Chisq")

}

#'  Print numeric with attributes
#'
#'  Prints generic named vectors without printing their attributes

#' @param object object of class "nwa" (numeric with attributes)
#'          or \code{\link[stats]{deviance}} is defined.
#' @param ... not implemented yet
#' @author Marcello Gallucci
#' @examples
#' data(manymodels)
#' model<-glm(ybin~x,family=binomial(),data=manymodels)
#' print(r2(model))
#' @export

print.nwa<-function(x,...) {
  a<-x
  attr(a,"df")<-NULL
  print(unclass(a))
}

