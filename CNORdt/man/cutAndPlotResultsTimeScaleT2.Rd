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
CNOlist, indexList, boolUpdates, divTime, lowerB=lowerB, upperB=upperB,
show = TRUE, plotPDF = FALSE, tag = NULL)
}
\arguments{
\item{model}{
The expanded model used as input for gaBinaryTimeScale.
}
\item{bStringT1}{
A bit string for T1 as output by gabinaryTimeScale (i.e. a vector of 1s and 0s).
}
\item{bStringT2}{
A bit string for T2 as output by gabinaryTimeScaleT2 (i.e. a vector of 1s and 0s).
}
\item{simList}{
A simlist corresponding to the model, as output by prep4sim.
}
\item{CNOlist}{
The CNOlist used in optimization.
}
\item{indexList}{
An indexList, produced by indexFinder run on the model and the CNOlist above.
}
\item{boolUpdates}{
The number of synchronous updates performed by the boolean simulator.
}
\item{divTime}{
The experimental time point that divides the early and late time phases, if using the 2-time phase version.
}
\item{lowerB}{
The lower bound for the optimized value of the scaling factor.
}
\item{upperB}{
The upper bound for the optimized value of the scaling factor.
}
\item{show}{
TOADD
}
\item{plotPDF}{
TOADD
}
\item{tag}{
TOADD
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
data(modelPB, package="CNORdt")

# pre-process model
model <- preprocessing(CNOlistPB, modelPB)
fields4Sim <- prep4sim(model=model)
initBstring <- rep(1, length(model$reacID))
indices = indexFinder(CNOlistPB, model)
# define boolUpdates and divTime
boolUpdates=c(10,20)
divTime=10

# optimize T1
opt1 <- gaBinaryTimeScale(CNOlist=CNOlistPB, model=model,
simList=fields4Sim, indexList=indices, initBstring=initBstring,
verbose=TRUE, boolUpdates=boolUpdates, divTime=divTime, maxTime=10,
lowerB=0.8, upperB=10)

# simulate the above model to get the starting point for second time phase
dataStartPoint = cutModel(model=model, simList=fields4Sim, bitString=opt1$bString)
simT1 = simulatorTimeScale(CNOlist=CNOlistPB, model=dataStartPoint[[1]], simList=dataStartPoint[[2]],
indexList=indices, boolUpdates=boolUpdates[1])

# optimize the 'second' network then stitch the results together
opt2 <- gaBinaryTimeScaleT2(CNOlist=CNOlistPB, model=model, simList=fields4Sim,
indexList=indices, bStringT1=opt1$bString, simResT1=simT1, verbose=TRUE,
boolUpdates=boolUpdates, divTime=divTime, lowerB=0.8, upperB=10)

# visualize total result
cutAndPlotResultsTimeScaleT2(model=model, bStringT1=opt1$bString, bStringT2=opt2$bString,
simList=fields4Sim, CNOlist=CNOlistPB, indexList=indices, boolUpdates=boolUpdates,
divTime=divTime, lowerB=0.8, upperB=10)
}
