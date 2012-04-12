\name{Binference}
\alias{Binference}

\title{
  Bayesian network inference
}

\description{
  This function uses data (CNOlist) to infer a Bayesian network using the catnet package.
}

\usage{
Binference(CNOlist, Model="AIC", tempCheckOrders=10,
            maxIter=100, filename="BAYESIAN")
}


\arguments{
  \item{CNOlist}{a CNOlist structure, as produced by \link{makeCNOlist}}
  \item{mode}{a character, optimization network selection criterion such as "AIC" and "BIC", to be used in \link{cnSearchSA}}
  \item{tempCheckOrders}{an integer, the number of iteration, orders to be searched, with constant temperature, to be used in \link{cnSearchSA}}
  \item{maxIter}{an integer, the total number of iterations, thus orders, to be processed, to be used in \link{cnSearchSA}}
  \item{filename}{name of the sif file saved, default BAYESIAN}
}
\details{
This function transforms the data in a format compatible with catnet package, infers the network using the Stochastic Network Search as implemented in catnet (see \link{cnSearchSA}), computes the consensus model of the models returned by \link{cnSearchSA} considering only links that have a frequency of appearence greater than 0.1 and returns the model in the sif format.
}

\value{
\item{sif}{the inferred data-driven network in sif format}
}

\author{
    F.Eduati
}

\seealso{
\link{MapDDN2Model}
}

\examples{

library(CellNOptR)
data(CNOlistDREAM,package="CellNOptR")
DDN<-Binference(CNOlistDREAM, tempCheckOrders=10, maxIter=100,
                filename="BAYESIAN")

    
}

