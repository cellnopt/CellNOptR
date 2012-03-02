\name{simulatorTimeScale}
\alias{simulatorTimeScale}
\title{
Discrete time simulation of a boolean model
}
\description{
Simulates multiple time points
}
\usage{
simulatorTimeScale(CNOlist, Model, SimList, indexList, boolUpdates)
}
\arguments{
\item{CNOlist}{
a CNOlist.
}
\item{Model}{
a Model that only contains the reactions to be evaluated.
}
\item{SimList}{
a SimList as created by prep4Sim, that has also already been cut to contain only the reactions to be evaluated.
}
\item{indexList}{
an indexList as created by indexFinder.
}
\item{boolUpdates}{
The number of update rounds the simulator should run for.
}
}
\details{
To simulate a particular model according to a user-defined bit string, use cutModel() and then simulate from the output.
}
\value{
A 3-dimensional array that gives the value of all species under each condition at each update (conditions, species, update).
}
\references{
}
\author{
A. MacNamara
}
\seealso{
gaBinaryTimeScale, getFitTimeScale
}
\examples{
# this computes the output of the full model, 
# which is normally not done on a stand alone basis,
# but if you have a model and would like to visualise 
# its output compared to your data, then this is what you should do.

library(CellNOptR)
library(CNORdt)
data(CNOlistToyFB, package="CNORdt")
data(ToyModelFB, package="CNORdt")

indexOrig <- indexFinder(CNOlistToyFB, ToyModelFB, verbose=TRUE)
fields4Sim <- prep4Sim(ToyModelFB)

simRes <- simulatorTimeScale(
  CNOlist=CNOlistToyFB,
  Model=ToyModelFB,
  SimList=fields4Sim,
  indexList=indexOrig,
  boolUpdates=20
)
}