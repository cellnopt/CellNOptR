\name{getFitTimeScale}
\alias{getFitTimeScale}
\title{
The optimization function that finds the scaling factor for the boolean simulation
}
\description{
This function is called from gaBinaryTimeScale. Using the model passed as input, it finds a scaling factor that minimizes the mean squared error between the data from the boolean simulation and the experimental data. A spline is fitted to the experimental data to allow this.
}
\usage{
getFitTimeScale(simList, CNOlist, model, indexList, sizeFac=0.0001,
NApenFac=1, boolUpdates, timeSplit="early", divTime=NULL, simResultsT1, lowerB=lowerB, upperB=upperB)
}
\arguments{
\item{simList}{
A list that contains additional fields for the simulator, as created by prep4sim applied to the model above.
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
\item{boolUpdates}{
The number of synchronous updates performed by the boolean simulator.
}
\item{timeSplit}{
If using the 2 steady-state version, this is the user-provided time point that signifies the split of early and late edges in the network.
}
\item{divTime}{
The experimental time point that divides the early and late steady states, if using the 2-steady state version.
}
\item{simResultsT1}{
The results of optimization at the first steady state (i.e. the whole list as output by gaBinaryTimeScale).
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
gaBinaryTimeScale, simulatorTimeScale
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

getFitRes <- getFitTimeScale(
  simList=fields4Sim,
  CNOlist=CNOlistPB,
  model=modelPB,
  indexList=indexOrig,
  boolUpdates=30,
  lowerB=0.8,
  upperB=10
)
}
