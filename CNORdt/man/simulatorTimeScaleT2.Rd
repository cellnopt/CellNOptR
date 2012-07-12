\name{simulatorTimeScaleT2}
\alias{simulatorTimeScaleT2}
\title{
Simulation of a boolean time course model for time 2
}
\description{
This function simulates a boolean model at time 2 where time 2 is the late time phase (late reactions) chosen by the user.
}
\usage{
simulatorTimeScaleT2(simResultsT1, CNOlist, model, simList, indexList, boolUpdates)
}
\arguments{
\item{simResultsT1}{
a matrix that is the output of simulatorT1 (i.e. one row per condition and one column per species).
}
\item{CNOlist}{
a CNOlist.
}
\item{model}{
a Model that only contains the reactions to be evaluated.}
\item{simList}{
a SimList as created by prep4sim, that has also already been cut to contain only the reactions to be evaluated.
}
\item{indexList}{
an indexList as created by indexFinder.
}
\item{boolUpdates}{
the number of synchronous updates performed by the boolean simulator.
}
}
\details{
Please see simulatorT2 for full details (simulatorTimeScaleT2 uses the same principles except it returns the full boolean simulation, not just the steady-state values).
}
\value{
A 3-dimensional array that gives the value of all species under each condition at each update (conditions, species, update).
}
\author{
A. MacNamara
}
\seealso{
gaBinaryTimeScaleT2, getFitTimeScale
}
\examples{
# this computes the output of the full model, 
# which is normally not done on a stand alone basis,
# but if you have a model and would like to visualize 
# its output compared to your data, then this is what you should do.

library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

indexOrig <- indexFinder(CNOlist, Model, verbose=TRUE)
fields4Sim <- prep4sim(Model)

simRes <- simulatorTimeScale(
  CNOlist=CNOlist,
  model=Model,
  simList=fields4Sim,
  indexList=indexOrig,
  boolUpdates=10
)

SimRes2 <- simulatorTimeScaleT2(
  simResultsT1=simRes[,,dim(simRes)[3]], 
  CNOlist=CNOlist, 
  model=Model, 
  simList=fields4Sim,
  indexList=indexOrig,
  boolUpdates=20)
}
