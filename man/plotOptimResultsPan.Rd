\name{plotOptimResultsPan}
\alias{plotOptimResultsPan}
\title{
Plot the data and simulated values
}
\description{
This function is the equivalent of CNOPlotFits, it plots the data and the simulated values,  along with an image plot that tells which cues were present.  The plots are coloured according to the fit between data and simulated data.
}
\usage{
plotOptimResultsPan(simResults, yInterpol=NULL, xCoords=NULL,
    CNOlist=CNOlist, nsplit=1, formalism=c("ss1","ss2","ssN","dt","ode"), pdf=FALSE,
    pdfFileName="", tPt=NULL)


}
\arguments{
  \item{simResults}{
a list with a field for each time point, each containing a matrix of dimensions (number of conditions) * (number of signals), with the first field being t0.  Typically produced by simulating a model and then extracting the columns that correspond to signals
}
  \item{yInterpol}{
todo
}
  \item{xCoords}{
todo
}

  \item{CNOlist}{
todo
}
  \item{nsplit}{
todo
}
  \item{formalism}{
todo
}
  \item{pdf}{
todo
}
  \item{pdfFileName}{
todo
}
  \item{tPt}{
todo
}
}
\details{
The colouring of the background is done as follows: the mean absolute difference
between observed and simulated values are computed, and colours are chosen based
on this value: red (1), yellow/orange (0.5),  white(0). This function is used inside \code{cutAndPlot}.
}
\value{
This function doesn't return anything, it just produces a plot in your graphics window.
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
A. MacNamara
}
\seealso{
cutAndPlotResultsT1, plotOptimResults
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)

#We will plot the fit of the full initial model compared to the data, without any optimisation
#This is normally not done on a stand alone basis, but if you have a model and would like to visualise its output compared to your data, then this is what you should do

#load and prepare data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

simList = prep4sim(ToyModel)
indexList = indexFinder(CNOlistToy, ToyModel)
#simulate model

#Sim<-simulatorTN(CNOlist=CNOlistToy,model=ToyModel, indexList=indexList, simList=simList)


## keep simList and indxList for back compatibility ?
#modelCut <- cutModel(model, bString)
#simListCut <- cutSimList(simList, bString)
#
    # t0
#Sim0 <- simulatorT0(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
#simRes0 <- as.matrix(Sim0[,indexList$signals])

# t1
#Sim <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
#simResults <- as.matrix(Sim[,indexList$signals])
#simResults <- list(t0=simRes0, t1=simRes)
#plot

#plotOptimResultsPan(
#	simResults=simResults,
#    CNOlist=CNOlistToy,
#    formalism="ss1",
#    tPt=1)
}
