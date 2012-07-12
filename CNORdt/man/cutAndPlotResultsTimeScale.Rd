\name{cutAndPlotResultsTimeScale}
\alias{cutAndPlotResultsTimeScale}
\title{
Plot the results of a time-course optimization
}
\description{
This function takes the optimized bit string, cuts the model according to the string and then finds the optimized scaling factor to pass on to plotOptimResultsTimeScale for visualization.
}
\usage{
cutAndPlotResultsTimeScale(model, bString, simList, CNOlist,
indexList, boolUpdates=boolUpdates, divTime=NULL, lowerB=lowerB, upperB=upperB)
}
\arguments{
\item{model}{
the expanded model used as input for gaBinaryTimeScale.
}
\item{bString}{
a bit string for T1 as output by gabinaryTimeScale (i.e. a vector of 1s and 0s).
}
\item{simList}{
a simlist corresponding to the model, as output by prep4sim.
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
\item{lowerB}{
The lower bound for the optimized value of the scaling factor.
}
\item{upperB}{
The upper bound for the optimized value of the scaling factor.
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
data(CNOlistPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

# pre-process model
indexOrig <- indexFinder(CNOlist=CNOlist, model=Model, verbose=TRUE)
indexNONC <- findNONC(model=Model, indexes=indexOrig, verbose=TRUE)
ModelCut <- cutNONC(model=Model, NONCindexes=indexNONC)
indexNONCcut <- indexFinder(CNOlist=CNOlist, model=ModelCut)
ModelCutCompress <- compressModel(model=ModelCut, indexes=indexNONCcut)
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, model=ModelCutCompress)
ModelCutCompressExpand <- expandGates(model=ModelCutCompress)
fields4Sim <- prep4sim(model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

# optimize
opt1 <- gaBinaryTimeScale(CNOlist=CNOlist, model=ModelCutCompressExpand,
simList=fields4Sim, indexList=indexNONCcutComp, initBstring=initBstring,
verbose=TRUE, boolUpdates=30, maxTime=100, lowerB=0.8, upperB=10)

#plot
cutAndPlotResultsTimeScale(
  model=ModelCutCompressExpand,
  bString=opt1$bString,
  simList=fields4Sim,
  CNOlist=CNOlist,
  indexList=indexNONCcutComp,
  boolUpdates=30,
  lowerB=0.8,
  upperB=10)
}
