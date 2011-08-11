\name{plotFit}
\alias{plotFit}
\title{
Plot the evolution of an optimisation
}
\description{
This function takes in the results of an optimisation by gaBinaryT1 and plots the evolution of best fit and average fit against generations.
}
\usage{
plotFit(OptRes)
}
\arguments{
  \item{OptRes}{
an object created by the optimisation engine (gabinaryT1)
}
}
\value{
This function doesn't return anything, it just produces a plot in your graphics window.
}
\author{
C. Terfve
}
\seealso{
gabinaryT1
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=FALSE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
indicesToyNCNOcut<-indexFinder(CNOlistToy,ToyNCNOcut)
ToyNCNOcutComp<-compressModel(ToyNCNOcut,indicesToyNCNOcut)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp)
ToyFields4Sim<-prep4Sim(ToyNCNOcutCompExp)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))
ToyT1opt<-gaBinaryT1(CNOlist=CNOlistToy,Model=ToyNCNOcutCompExp,SimList=ToyFields4Sim,indexList=indicesToyNCNOcutComp,initBstring=initBstring,maxGens=3,verbose=FALSE)
plotFit(OptRes=ToyT1opt)
}
