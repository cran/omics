\name{mlmer}
\alias{mlmer}
\title{Multiple Linear Mixed-Effects Models}
\description{
    Fits multiple linear mixed-effects models.
}
\usage{
    mlmer(formula, data, vars, lrt=TRUE, save.residuals=FALSE, save.ranks=TRUE)
}
\arguments{
    \item{formula}{an object of class \code{formula} (or one that can be coerced
                   to that class): a symbolic description of the models to be
                   fitted.}
    \item{data}{an optional data frame, list or environment (or object coercible
                by \code{as.data.frame} to a data frame) containing the
                variables in the model.}
    \item{vars}{a character vector of variables of interest.}
    \item{lrt}{whether \eqn{p}-values should be computed using likelihood-ratio
               tests.}
    \item{save.residuals}{whether model residuals should be stored and
                          returned.}
    \item{save.ranks}{whether random effect ranks should be stored and
                      returned.}
}
\value{
    If \code{save.ranks} is \code{TRUE} (default), a list with elements:
    \item{coefs}{an array (or data frame, if \code{vars} contains a single
                 element) of regression coefficient estimates (\code{coef}),
                 standard errors (\code{coef.se}), and corresponding
                 \eqn{p}-values (\code{pval}).}
    \item{ranef.ranks}{a list of rank matrices (one per random effect) that can
                       be visualized using \code{\link{ranks.heatmap}}.}

    If \code{save.residuals} is additionally \code{TRUE}, the list also
    includes:
    \item{residuals}{a matrix of model residuals.}

    If both \code{save.ranks} and \code{save.residuals} are \code{FALSE}, only
    \code{coefs} is returned (similarly to \code{\link{mlm}}).
}
\seealso{
    \code{\link{lmer}} and \code{\link{mlm}}
}
