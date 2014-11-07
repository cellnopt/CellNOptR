\name{CNORfeeder-package}
\alias{CNORfeeder-package}
\alias{CNORfeeder}
\docType{package}
\title{
R package to integrate literature-constrained and data-driven methods to infer signalling networks from perturbation experiments
}
\description{
CNORfeeder permits to extend a network derived from literature with links derived strictly from the data via various inference methods using information on physical interactions of proteins to guide and validate the integration of links. The package is designed to be integrated with CellNOptR.
}
\details{
\tabular{ll}{
Package: \tab CNORfeeder\cr
Type: \tab Package\cr
Version: \tab 1.0.0.\cr
Date: \tab 2012-11-22\cr
License: \tab GPLv2\cr
LazyLoad: \tab yes\cr
}
}
\author{
F. Eduati
Maintainer: F. Eduati <eduati@ebi.ac.uk>
}
\references{
F. Eduati, J. De Las Rivas, B. Di Camillo, G. Toffolo, J. Saez-Rodriguez. Integrating literature-constrained and data-driven inference of signalling networks. Bioinformatics, 28(18):2311-2317, 2012.
}

\keyword{ package }
\examples{
library(CNORfeeder)

# this is an example of the main steps of the integrated CellNOptR - CNORfeeder pipeline

# load the data already formatted as CNOlist
data(CNOlistDREAM,package="CellNOptR")
# load the model (PKN) already in the CNO format
data(DreamModel,package="CellNOptR")
# see CellNOptR documentation to import other data/PKNs)

# A. INFERENCE - CNORfeeder
# FEED inference: codified in Boolean Tables
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))

# B. COMPRESSION - CellNOptR
# preprocessing step
model<-preprocessing(data=CNOlistDREAM, model=DreamModel)

# C. INTEGRATION - CNORfeeder
# integration with the compressed model
modelIntegr <- mapBTables2model(BTable=BTable,model=model,allInter=TRUE)
# see example in ?MapDDN2Model to use other reverse-engineering methods

# D. WEGHTING - CNORfeeder
# integrated links are weighted more according to the integratin factor integrFac
modelIntegrWeight <- weighting(modelIntegr=modelIntegr, PKNmodel=DreamModel,
                              CNOlist=CNOlistDREAM, integrFac=10)

# E. TRAINING - CellNOptR
initBstring<-rep(1,length(modelIntegr$reacID))
# training to data using genetic algorithm (run longer to obtain better results)
DreamT1opt<-gaBinaryT1W(
CNOlist=CNOlistDREAM,
	model=modelIntegrWeight,
	initBstring=initBstring,
	maxGens=2,
	popSize=5,
	verbose=FALSE)
}