\name{simulatorDT}
\alias{simulatorDT}
\title{
Discrete time simulation of a boolean model
}
\description{
Simulates multiple time points within C (for speed).
}
\usage{
simulatorDT(CNOlist, model, simList, indexList, boolUpdates)
}
\arguments{
\item{CNOlist}{
A CNOlist.
}
\item{model}{
A model that only contains the reactions to be evaluated.
}
\item{simList}{
A simList as created by prep4sim, that has also already been cut to contain only the reactions to be evaluated.
}
\item{indexList}{
An indexList as created by indexFinder.
}
\item{boolUpdates}{
The number of update rounds the simulator should run for.
}
\item{prevSim}{
The results from \code{simulatorT0} can be used here as initial conditions.
}
}
\details{

}
\value{
A 3-dimensional array that gives the value of all species under each condition at each update (conditions, species, update).
}
\author{
A. MacNamara
}
\seealso{
gaBinaryDT, getFitDT
}
\examples{
# this computes the output of the full model, 
# which is normally not done on a stand alone basis,
# but if you have a model and would like to visualise 
# its output compared to your data, then this is what you should do.

library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

indexOrig <- indexFinder(CNOlistPB, modelPB, verbose=TRUE)
fields4Sim <- prep4sim(modelPB)

simResults <- simulatorDT(
  CNOlist=CNOlistPB,
  model=modelPB,
  simList=fields4Sim,
  indices=indexOrig,
  boolUpdates=10
)
simResults = convert2array(simResults, dim(CNOlistPB$valueSignals[[1]])[1],
length(modelPB$namesSpecies), boolUpdates)
}
