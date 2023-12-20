#'  LRT Power
#'
#' computes power parameters based on the R-squared or eta-squared for
#' generalized linear model with categorical dependent variables

#' @param es The expected r-square or expected eta-square
#' @param df degrees of freedom of the model or of the effect
#' @param prob expected proportions of cases in the dependent variables levels. For dichotomous
#'            dependent variables a number between 0 and 1 (excluded). For categorical variables with more
#'            than two levels, a numeric vector of length equal to the number of levels and whose
#'            values sum to one.
#' @param N  Number of observations. Omit if required sample size N is required
#' @param sig.level Significance level (Type I error probability). Set to NULL if required minimal alpha is required
#' @param power Power of test (1 minus Type II error probability). Omit if estimated power is required.

#'
#' @return Object of class '"power.htest"' with input and computed parameters
#' @author Marcello Gallucci
#' @examples
#' @rdname power
#' @export

power.lrt <- function(es,prob,df=NULL, N=NULL,sig.level=.05, power=NULL) {

  if (is.null(prob))
    stop("Dependent variable levels proportion should be defined")

  if (sum(sapply(list(es, N, df, power, sig.level), is.null)) !=  1)
      stop("exactly one of es, N, df, power, and sig.level must be NULL")
  if (length(prob)==1) prob<-c(prob,1-prob)
  if (sum(prob)!=1)
      stop("Dependent variable levels proportion should sum to 1")
  D0<--2*sum(prob*log(prob))
  lambda<-sqrt(es*D0)
  if (length(lambda)==0) lambda<-NULL
  res<-pwr::pwr.chisq.test(lambda,df=df,power = power,sig.level = sig.level,N=N)
  if (is.null(es))
      res$w<-res$w^2/D0
  else
      res$w<-es

  names(res)[1]<-"es"
  if ((length(prob)-1)>df)
    message("DF are less than the number of dummies representing the dependent variable. Please be sure that the df are correct.")
  res

}


