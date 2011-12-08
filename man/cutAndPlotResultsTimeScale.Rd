\name{cutAndPlotResultsTimeScale}
\alias{cutAndPlotResultsTimeScale}
\title{
Plot the results of a time-course optimisation
}
\description{
This function takes the optimized bit string, cuts the model according to the string and then finds the optimized scaling factor to pass on to plotOptimResultsTimeScale for visualization.
}
\usage{
cutAndPlotResultsTimeScale(Model, bString, SimList, CNOlist,
indexList, boolUpdates=boolUpdates, divTime=NULL)
}
\arguments{
\item{Model}{
a model (the full one that was used for optimisation)
}
\item{bString}{
a bit string for T1 as output by gabinaryT1 (i.e. a vector of 1s and 0s)
}
\item{SimList}{
a simlist corresponding to the model, as output by prep4Sim
}
\item{CNOlist}{
the CNOlist used in optimization
}
\item{indexList}{
an indexList, produced by indexFinder run on the model and the CNOlist above
}
\item{boolUpdates}{
the number of synchronous updates performed by the boolean simulator.
}
\item{divTime}{
the experimental time point that divides the early and late steady states, if using the 2-steady state version.
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
data(CNOlistToyFB, package="CNORdiscreteTime")
data(ToyModelFB, package="CNORdiscreteTime")

# pre-process model

indexOrig <- indexFinder(CNOlist=CNOlistToyFB, Model=ToyModelFB, verbose=TRUE)
indexNONC <- findNONC(Model=ToyModelFB, indexes=indexOrig, verbose=TRUE)
ModelCut <- cutNONC(Model=ToyModelFB, NONCindexes=indexNONC)
indexNONCcut <- indexFinder(CNOlist=CNOlistToyFB, Model=ModelCut)
ModelCutCompress <- compressModel(Model=ModelCut, indexes=indexNONCcut)
indexNONCcutComp <- indexFinder(CNOlist=CNOlistToyFB, Model=ModelCutCompress)
ModelCutCompressExpand <- expandGates(Model=ModelCutCompress)
fields4Sim <- prep4Sim(Model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

# optimize

boolUpdates = 20
optFB <- gaBinaryTimeScale(
  CNOlist=CNOlistToyFB,
  Model=ModelCutCompressExpand,
  SimList=fields4Sim,
  indexList=indexNONCcutComp,
  initBstring=initBstring,
  verbose=TRUE,
  boolUpdates=boolUpdates,
  divTime=NULL)

cutAndPlotResultsTimeScale(
  Model=ModelCutCompressExpand,
  bString=optFB$bString,
  SimList=fields4Sim,
  CNOlist=CNOlistToyFB,
  indexList=indexNONCcutComp,
  boolUpdates=boolUpdates)
}

