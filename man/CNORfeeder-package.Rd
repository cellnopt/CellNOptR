\name{CNORfeeder-package}
\alias{CNORfeeder-package}
\alias{CNORfeeder}
\docType{package}
\title{
R version of CNORfeeder, boolean features
}
\description{
This package does optimisation of boolean logic networks of signalling pathways based on a previous knowledge network and a set of data collected upon perturbation of some of the nodes in the network.
}
\details{
\tabular{ll}{
Package: \tab CNORfeeder\cr
Type: \tab Package\cr
Version: \tab 1.5.6\cr
Date: \tab 2012-11-22\cr
License: \tab GPLv2\cr
LazyLoad: \tab yes\cr
}
}
\author{
T. Cokelaer, A. MacNamara, F. Eduati, S. Schrier, C. Terfve

Maintainer: T. Cokelaer <cokelaer@ebi.ac.uk>
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}

\keyword{ package }
\examples{
library(CNORfeeder)

# quick 1 time point optimisation of a Prior Knowledge Network to MIDAS data. 
data(CNOlistToy,package="CNORfeeder")
data(ToyModel,package="CNORfeeder")
model = preprocessing(CNOlistToy, ToyModel)
results = gaBinaryT1(CNOlistToy, model, verbose=FALSE)
plotFit(results)
cutAndPlot(CNOlistToy, model, list(results$bString))

# Same as above and HTML report
CNORwrap(name="Toy", 
    namesData=list(CNOlist="ToyData",model="ToyModel"),
    data=CNOlistToy, model=ToyModel)
}

