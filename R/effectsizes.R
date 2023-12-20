
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

r2.default<-function(object,...) {


  d1<-stats::deviance(object)
  if (is.null(d1)) {
    d1<- -2*stats::logLik(object)
  }
  model0<-stats::update(object,.~1)
  df<-length(stats::coef(object))-length(stats::coef(model0))
  d0<-stats::deviance(model0)
  if (is.null(d0))
    d0<- -2*stats::logLik(model0)

  r2<-1-d1/d0
  r2adj<-1-(d1+df)/d0
  obj<-c(r2=r2,r2adj=r2adj)
  attr(obj,"df")<-df
  obj
}

#' @rdname r2
#' @export

r2.lm<-function(object,...) {

  ss<-summary(object)
  c(r2=  ss$r.squared,r2adj=  ss$adj.r.squared)
}

#' @rdname r2
#' @export

r2.glm<-function(object,...) {

  r2.default(object)
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

  model0<-stats::update(object ,.~.)
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
  model0<-stats::update(object ,.~1)
  w<-which(rownames(a) %in% c("(Intercept)","Residuals"))
  res<-a$`Sum Sq`/(model0$df.residual*stats::sigma(model0)^2)
  res<-matrix(res[-w],ncol=1)
  rownames(res)<-rownames(a)[-w]
  colnames(res)<-"Eta_squared"

  eps<-(a$`Sum Sq`-a$`Sum Sq`[length(a$`Sum Sq`)]/object$df.residual)/(model0$df.residual*stats::sigma(model0)^2)
  eps<-matrix(eps[-w],ncol=1)
  rownames(eps)<-rownames(a)[-w]
  colnames(eps)<-"Epsilon_squared"
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

  model0<-stats::update(object ,.~.)
  dev0<-stats::deviance(model0)
  if (is.null(dev0))
    dev0<- as.numeric(-2*stats::logLik(model0))
  a<-car::Anova(object,type=3,test=test)
  df<-a$Df
  print(a)
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
