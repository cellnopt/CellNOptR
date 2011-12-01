\name{plotOptimResultsTimeScale}
\alias{plotOptimResultsTimeScale}
\title{
Plot the data and simulated values for time course data
}
\description{
This function is the equivalent of CNOPlotFits, it plots the data and the simulated values,  along with an image plot that tells which cues were present. The plots are colored according to the fit between data and simulated data.
}
\usage{
plotOptimResults(SimResults = SimResults, yInterpol = yInterpol, xCoords = xCoords,
times = times, namesCues = namesCues, namesSignals = namesSignals, valueCues = valueCues)
}
\arguments{
\item{SimResults}{
a list with a field for each time point, each containing a matrix of dimensions (number of conditions) * (number of signals), with the first field being t=0. Typically produced by simulating a model and then extracting the columns that correspond to signals
}
\item{yInterpol}{
the interpolated experimental data returned by \code{getFitTimeScale}
}
\item{xCoords}{
the x-coordinates of the data and simulation after the application of the optimized scaling factor
}
\item{namesCues}{
a vector of names, typically \code{CNOlist$namesCues}
}
\item{namesSignals}{
a vector of names, typically \code{CNOlist$namesSignals}
}
\item{valueCues}{
a matrix of dimensions (number of conditions) * (number of cues), typically \code{CNOlist$valueCues}
}
}
\details{
A single coloring denotes the goodness of fit between simulation and experiment, this needs to be looked at again.
}
\value{
This function doesn't return anything, it just produces a plot in your graphics window.
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
A. MacNamara
}
\seealso{
cutAndPlotResultsTimeScale
}
\examples{
# this function is usually contained within cutAndPlotResultsTimeScale
# but the output can be viewed as follows:

library(CellNOptR)
library(CNORdiscreteTime)
data(CNOlistToyFB, package="CNORdiscreteTime")
data(ToyModelFB, package="CNORdiscreteTime")
boolUpdates=20

indexOrig <- indexFinder(CNOlistToyFB, ToyModelFB, verbose=TRUE)
fields4Sim <- prep4Sim(ToyModelFB)

SimRes <- simulatorTimeScale(
  CNOlist=CNOlistToyFB,
  Model=ToyModelFB,
  SimList=fields4Sim,
  indexList=indexOrig,
  boolUpdates=boolUpdates
)

SimRes <- SimRes[,indexOrig$signals,]

getFitRes <- getFitTimeScale(
 SimList=fields4Sim,
 CNOlist=CNOlistToyFB,
 Model=ToyModelFB,
 indexList=indexOrig,
 boolUpdates=boolUpdates
)

plotOptimResultsTimeScale(
  SimResults=SimRes,
  yInterpol=getFitRes$yInter,
  xCoords=getFitRes$xCoords,
  expResults=CNOlistToyFB$valueSignals,
  times=CNOlistToyFB$timeSignals,
  namesCues=CNOlistToyFB$namesCues,
  namesSignals=CNOlistToyFB$namesSignals,
  valueCues=CNOlistToyFB$valueCues)
}
