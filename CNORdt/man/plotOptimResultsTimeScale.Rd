\name{plotOptimResultsTimeScale}
\alias{plotOptimResultsTimeScale}
\title{
Plot the data and simulated values for time course data
}
\description{
This function is the equivalent of CNOPlotFits, it plots the data and the simulated values,  along with an image plot that tells which cues were present. The plots are colored according to the fit between data and simulated data.
}
\usage{
plotOptimResultsTimeScale(SimResults, yInterpol, xCoords, CNOlist, nsplit=1)
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
\item{CNOlist}{
a CNOlist on which the score is based (based on all valueSignals).
}
\item{nsplit}{
a value that splits the plot 'nsplit' times.
}
}
\details{
The fit is colored using a heatmap from white (perfect fit) to red.
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
data(CNOlistToyPB, package="CNORdiscreteTime")
data(ToyModelPB, package="CNORdiscreteTime")
boolUpdates=20

indexOrig <- indexFinder(CNOlist, Model, verbose=TRUE)
fields4Sim <- prep4Sim(Model)

SimRes <- simulatorTimeScale(
   CNOlist=CNOlist,
   Model=Model,
   SimList=fields4Sim,
   indexList=indexOrig,
   boolUpdates=boolUpdates
)

SimRes <- SimRes[,indexOrig$signals,]

getFitRes <- getFitTimeScale(
  SimList=fields4Sim,
  CNOlist=CNOlist,
  Model=Model,
  indexList=indexOrig,
  boolUpdates=boolUpdates,
  lowerB=0.8,
  upperB=10
)

plotOptimResultsTimeScale(
   SimResults=SimRes,
   yInterpol=getFitRes$yInter,
   xCoords=getFitRes$xCoords,
   CNOlist=CNOlist,
   nsplit=1
)
}
