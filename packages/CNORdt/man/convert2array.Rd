\name{convert2array}
\alias{convert2array}
\title{
Converts the output of simulatorDT to a 3D array for other functions
}
\description{
As above.
}
\usage{
convert2array(x, nRow, nCol, nBool)
}
\arguments{
\item{x}{
simResults as returned by simulatorDT
}
\item{nRow}{
Number of rows of output.
}
\item{nCol}{
Number of columns of output.
}
\item{nBool}{
The third dimension, which corresponds to the number of simulation updates.
}
}
\value{
This function returns a 3D array of the simulation results.
}
\author{
A. MacNamara
}
\examples{
library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

indexOrig <- indexFinder(CNOlistPB, modelPB, verbose=TRUE)
fields4Sim <- prep4sim(modelPB)
boolUpdates=10

simResults <- simulatorDT(
  CNOlist=CNOlistPB,
  model=modelPB,
  simList=fields4Sim,
  indices=indexOrig,
  boolUpdates=boolUpdates
)
simResults = convert2array(simResults, dim(CNOlistPB$valueSignals[[1]])[1],
length(modelPB$namesSpecies), boolUpdates)
}