\name{cutAndPlotResultsTimeScaleT2}
\alias{cutAndPlotResultsTimeScaleT2}
\title{
Plot the results of a 2 time-phase time-course optimization
}
\description{
This function takes the optimized bit string of the early and late time phases, cuts the model according to the strings and then finds the optimized scaling factor to pass on to plotOptimResultsTimeScale for visualization.
}
\usage{
cutAndPlotResultsTimeScaleT2(model, bStringT1, bStringT2, simList,
CNOlist, indexList, boolUpdates, divTime, lowerB=lowerB, upperB=upperB)
}
\arguments{
\item{model}{
the expanded model used as input for gaBinaryTimeScale.
}
\item{bStringT1}{
a bit string for T1 as output by gabinaryTimeScale (i.e. a vector of 1s and 0s).
}
\item{bStringT2}{
a bit string for T2 as output by gabinaryTimeScaleT2 (i.e. a vector of 1s and 0s).
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
library(CellNOptR)
library(CNORdt)
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

# define boolUpdates and divTime
boolUpdates=c(10,20)
divTime=10

# optimize T1
opt1 <- gaBinaryTimeScale(CNOlist=CNOlist, model=ModelCutCompressExpand,
simList=fields4Sim, indexList=indexNONCcutComp, initBstring=initBstring,
verbose=TRUE, boolUpdates=boolUpdates, divTime=divTime, maxTime=100,
lowerB=0.8, upperB=10)

# simulate the above model to get the starting point for second time phase
dataStartPoint = cutModel(model=ModelCutCompressExpand, simList=fields4Sim, bitString=opt1$bString)
simT1 = simulatorTimeScale(CNOlist=CNOlist, model=dataStartPoint[[1]], simList=dataStartPoint[[2]],
indexList=indexNONCcutComp, boolUpdates=boolUpdates[1])

# optimize the 'second' network then stitch the results together
opt2 <- gaBinaryTimeScaleT2(CNOlist=CNOlist, model=ModelCutCompressExpand,
simList=fields4Sim, indexList=indexNONCcutComp,
bStringT1=opt1$bString, simResT1=simT1, verbose=TRUE, boolUpdates=boolUpdates, divTime=divTime, lowerB=0.8, upperB=10)

# visualize total result
cutAndPlotResultsTimeScaleT2(model=ModelCutCompressExpand, bStringT1=opt1$bString, bStringT2=opt2$bString,
simList=fields4Sim, CNOlist=CNOlist, indexList=indexNONCcutComp, boolUpdates=boolUpdates,
divTime=divTime, lowerB=0.8, upperB=10)
}
