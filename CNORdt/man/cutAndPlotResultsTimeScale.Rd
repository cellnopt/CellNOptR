\name{cutAndPlotResultsTimeScale}
\alias{cutAndPlotResultsTimeScale}
\title{
Plot the results of a time-course optimization
}
\description{
This function takes the optimized bit string, cuts the model according to the string and then finds the optimized scaling factor to pass on to plotOptimResultsTimeScale for visualization.
}
\usage{
cutAndPlotResultsTimeScale(Model, bString, SimList, CNOlist,
indexList, boolUpdates=boolUpdates, divTime=NULL, lowerB=lowerB, upperB=upperB)
}
\arguments{
\item{Model}{
the expanded model used as input for gaBinaryTimeScale.
}
\item{bString}{
a bit string for T1 as output by gabinaryTimeScale (i.e. a vector of 1s and 0s).
}
\item{SimList}{
a simlist corresponding to the model, as output by prep4Sim.
}
\item{CNOlist}{
the CNOlist used in optimization.
}
\item{indexList}{
an indexList, produced by indexFinder run on the model and the CNOlist above.
}
\item{boolUpdates}{
the number of synchronous updates performed by the boolean simulator.
}
\item{divTime}{
the experimental time point that divides the early and late time phases, if using the 2-time phase version.
}
}

\value{
This function doesn't return anything, it only plots the graph in your graphic window.
}
\author{
A. MacNamara
}
\seealso{
plotOptimResultsTimeScale
}
\examples{
data(CNOlistToyPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

# pre-process model
indexOrig <- indexFinder(CNOlist=CNOlist, Model=Model, verbose=T)
indexNONC <- findNONC(Model=Model, indexes=indexOrig, verbose=T)
ModelCut <- cutNONC(Model=Model, NONCindexes=indexNONC)
indexNONCcut <- indexFinder(CNOlist=CNOlist, Model=ModelCut)
ModelCutCompress <- compressModel(Model=ModelCut, indexes=indexNONCcut)
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, Model=ModelCutCompress)
ModelCutCompressExpand <- expandGates(Model=ModelCutCompress)
fields4Sim <- prep4Sim(Model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

# optimize
opt1 <- gaBinaryTimeScale(CNOlist=CNOlist, Model=ModelCutCompressExpand,
SimList=fields4Sim, indexList=indexNONCcutComp, initBstring=initBstring,
verbose=TRUE, boolUpdates=30, MaxTime=100, lowerB=0.8, upperB=10)

#plot
cutAndPlotResultsTimeScale(
  Model=ModelCutCompressExpand,
  bString=opt1$bString,
  SimList=fields4Sim,
  CNOlist=CNOlist,
  indexList=indexNONCcutComp,
  boolUpdates=boolUpdates,
  lowerB=0.8,
  upperB=10)
}
