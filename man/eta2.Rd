% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effectsizes.R
\name{eta2}
\alias{eta2}
\alias{eta2.default}
\alias{eta2.lm}
\alias{eta2.glm}
\alias{eta2.clm}
\alias{eta2_partial.default}
\title{Eta-squared and Epsilon-squared}
\usage{
eta2(object, ...)

\method{eta2}{default}(object, ...)

\method{eta2}{lm}(object, ...)

\method{eta2}{glm}(object, ...)

\method{eta2}{clm}(object, ...)

\method{eta2_partial}{default}(object, ...)
}
\arguments{
\item{object}{object of class "glm", or an object for which the function \code{\link[stats]{drop1}}
or \code{\link[stats]{deviance}} is defined.}

\item{...}{not implemented yet}
}
\value{
an anova table with R
}
\description{
Computes the eta-squared and epsilon-squared indices for several generalized linear models
}
\examples{
data(manymodels)
model<-glm(ybin~x,family=binomial(),data=manymodels)
eta2(model)
}
\author{
Marcello Gallucci
}
