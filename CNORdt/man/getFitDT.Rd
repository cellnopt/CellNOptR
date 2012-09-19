\name{getFitTimeScale}
\alias{getFitTimeScale}
\title{
The optimization function that finds the scaling factor for the boolean simulation
}
\description{
This function is called from gaBinaryDT. Using the model passed as input, it finds a scaling factor that minimizes the mean squared error between the data from the boolean simulation and the experimental data. A spline is fitted to the experimental data to allow this.
}
\usage{
getFitDT(simResults, CNOlist, model, indexList, sizeFac=0.0001,
NApenFac=1, nInTot, boolUpdates, lowerB=lowerB, upperB=upperB)
}
\arguments{
\item{simResults}{
The simulation results as output from simulatorDT
}
\item{CNOlist}{
A CNOlist on which the score is based (based on all valueSignals).
}
\item{model}{
A model list.
}
\item{indexList}{
A list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on the model and CNOlist above.
}
\item{sizeFac}{
The scaling factor for the size term in the objective function, default to 0.0001.
}
\item{NApenFac}{
The scaling factor for the NA term in the objective function, default to 1.
}
\item{nInTot}{
The number of inputs in the model prior to cutting, used to normalise the size penalty.
}
\item{boolUpdates}{
The number of synchronous updates performed by the boolean simulator.
}
\item{lowerB}{
The lower bound for the optimized value of the scaling factor.
}
\item{upperB}{
The upper bound for the optimized value of the scaling factor.
}
}
\details{
The function optim() is used to find the optimal scaling factor.
}
\value{
This function returns a list with elements:
\item{score}{The mean squared error between simulation and experiment with NA and model size penalties.}
\item{estimate}{The scaling factor used to compare boolean and real data.}
\item{xCoords}{The x-axis coordinates after multiplication with the scaling factor.}
\item{yInter}{The interpolated values of the experimental data.}
\item{yBool}{The boolean simulation results at each time point.}
}
\author{
A. MacNamara
}
\seealso{
gaBinaryDT, simulatorDT
}
\examples{
# this function is usually contained within gaBinaryTimeScale
# but the output can be viewed as follows:

library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

# pre-processing
indexOrig <- indexFinder(CNOlist=CNOlistPB, model=modelPB, verbose=TRUE)
fields4Sim <- prep4sim(model=modelPB)

boolUpdates = 10
simResults <- simulatorDT(
  CNOlist=CNOlistPB,
  model=modelPB,
  simList=fields4Sim,
  indices=indexOrig,
  boolUpdates=boolUpdates
)
simResults = convert2array(simResults, dim(CNOlistPB$valueSignals[[1]])[1],
length(modelPB$namesSpecies), boolUpdates)

optimRes <- getFitDT(
  simResults=simResults,
  CNOlist=CNOlistPB,
  model=modelPB,
  indexList=indexOrig,
  boolUpdates=boolUpdates,
  lowerB=0.8,
  upperB=10,
  nInTot=length(which(modelPB$interMat == -1))
)
}
