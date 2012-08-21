\name{plotFit}
\alias{plotFit}
\title{
Plot the evolution of an optimisation
}
\description{
This function takes in the results of an optimisation by gaBinaryT1 and plots the evolution of best fit and average fit against generations.
}
\usage{
plotFit(optRes, filename = NULL)
}
\arguments{
  \item{optRes}{
an object created by the optimisation engine (gabinaryT1)
}
  \item{filename}{
NULL or string: if provided, the plot is save in PDF format in the filename.
}
}
\value{
This function doesn't return anything, it just produces a plot in your graphics window.
}
\author{
C. Terfve
}
\seealso{
gaBinaryT1
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#process the model

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyComp<-compressModel(ToyModel,indicesToy)
indicesToyComp<-indexFinder(CNOlistToy,ToyComp)
ToyCompExp<-expandGates(ToyComp)

#optimise

ToyFields4Sim<-prep4sim(ToyCompExp)
initBstring<-rep(1,length(ToyCompExp$reacID))
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	model=ToyCompExp,
	initBstring=initBstring,
	maxGens=2,
	popSize=5,
	verbose=FALSE)

plotFit(optRes=ToyT1opt)
}
