\name{getFitTimeScale}
\alias{getFitTimeScale}
\title{
The optimization function that finds the scaling factor for the boolean simulation
}
\description{
This function is called from gaBinaryTimeScale. Using the model passed as input, it finds a scaling factor that minimizes the mean squared error between the data from the boolean simulation and the experimental data. A spline is fitted to the experimental data to allow this.
}
\usage{
getFitTimeScale(SimList, CNOlist, Model, indexList, sizeFac=0.0001,
NAPenFac=1, boolUpdates, timeSplit="early", divTime=NULL, SimResultsT1, lowerB=lowerB, upperB=upperB)
}
\arguments{
\item{SimList}{
a list that contains additional fields for the simulator, as created by prep4Sim applied to the model above.
}
\item{CNOlist}{
a CNOlist on which the score is based (based on all valueSignals).
}
\item{Model}{
a Model list.
}
\item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on the model and CNOlist above.
}
\item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001.
}
\item{NAPenFac}{
the scaling factor for the NA term in the objective function, default to 1.
}
\item{boolUpdates}{
the number of synchronous updates performed by the boolean simulator.
}
\item{timeSplit}{
If using the 2 steady-state version, this is the user-provided time point that signifies the split of early and late edges in the network.
}
\item{divTime}{
the experimental time point that divides the early and late steady states, if using the 2-steady state version.
}
\item{SimResultsT1}{
the results of optimization at the first steady state (i.e. the whole list as output by gaBinaryTimeScale).
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
\item{score}{the mean squared error between simulation and experiment with NA and model size penalties.}
\item{estimate}{the scaling factor used to compare boolean and real data.}
\item{xCoords}{the x-axis coordinates after multiplication with the scaling factor.}
\item{yInter}{the interpolated values of the experimental data.}
\item{yBool}{the boolean simulation results at each time point.}
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
data(CNOlistToyPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

# pre-processing
indexOrig <- indexFinder(CNOlist=CNOlist, Model=Model, verbose=T)
fields4Sim <- prep4Sim(Model=Model)

getFitRes <- getFitTimeScale(
  SimList=fields4Sim,
  CNOlist=CNOlist,
  Model=Model,
  indexList=indexOrig,
  boolUpdates=30,
  lowerB=0.8,
  upperB=10
)
}
