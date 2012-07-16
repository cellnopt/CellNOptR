\name{CellNOptR-package}
\alias{CellNOptR-package}
\alias{CellNOptR}
\docType{package}
\title{
R version of CellNOptR, boolean features
}
\description{
This package does optimisation of boolean logic networks of signalling pathways based on a previous knowledge network and a set of data collected upon perturbation of some of the nodes in the network.
}
\details{
\tabular{ll}{
Package: \tab CellNOptR\cr
Type: \tab Package\cr
Version: \tab 1.3.23\cr
Date: \tab 2012-07-17\cr
License: \tab GPLv2\cr
LazyLoad: \tab yes\cr
}
}
\author{
C.Terfve, T.Cokelaer, A.MacNamara, F. Eduati


Maintainer: T.Cokelaer <cokelaer@ebi.ac.uk>
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}

\keyword{ package }
\examples{
library(CellNOptR)

#version with paramslist

tmpdir<-tempdir()
setwd(tmpdir)

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

pList = defaultParameters()
pList$data = CNOlistToy
pList$model = ToyModel
pList$maxGens = 5
pList$popSize = 5

CNORwrap(paramsList=pList, name="Toy", 
    namesData=list(CNOlist="ToyData", model="ToyModel"))

\dontrun{
#version with default parameters


CNORwrap(name="Toy", 
    namesData=list(CNOlist="ToyData",model="ToyModel"),
    data=CNOlistToy, model=ToyModel)
}
}
