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
indexList, boolUpdates=boolUpdates, divTime=NULL, lowerB=lowerB,
upperB=upperB, show = TRUE, plotPDF = FALSE, tag = NULL)
}
\arguments{
\item{model}{
The expanded model used as input for gaBinaryTimeScale.
}
\item{bString}{
A bit string for T1 as output by gabinaryTimeScale (i.e. a vector of 1s and 0s).
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
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

# pre-process model
model <- preprocessing(CNOlistPB, modelPB)
indices = indexFinder(CNOlistPB, model)
fields4Sim <- prep4sim(model=model)
initBstring <- rep(1, length(model$reacID))

# optimize
opt1 <- gaBinaryTimeScale(CNOlist=CNOlistPB, model=model,
simList=fields4Sim, indexList=indices, initBstring=initBstring,
verbose=TRUE, boolUpdates=30, maxTime=100, lowerB=0.8, upperB=10)

# plot
cutAndPlotResultsTimeScale(
  model=model,
  bString=opt1$bString,
  simList=fields4Sim,
  CNOlist=CNOlistPB,
  indexList=indices,
  boolUpdates=30,
  lowerB=0.8,
  upperB=10)
}
