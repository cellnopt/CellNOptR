\name{cutAndPlotResultsDT}
\alias{cutAndPlotResultsDT}
\title{
Plot the results of a time-course optimization
}
\description{
This function takes the optimized bit string, cuts the model according to the string and then finds the optimized scaling factor to pass on to plotOptimResultsPan for visualization.
}
\usage{
cutAndPlotResultsDT(model, bString, simList=NULL, CNOlist, indexList=NULL,
plotPDF=FALSE, tag=NULL, plotParams=list(maxrow=10), boolUpdates=boolUpdates,
lowerB=lowerB, upperB=upperB, sizeFac = 1e-04, NAFac = 1)
}
\arguments{
\item{model}{
The expanded model used as input for gaBinaryDT.
}
\item{bString}{
A bit string as output by gabinaryDT (i.e. a vector of 1s and 0s).
}
\item{simList}{
A simlist corresponding to the model, as output by prep4sim.
}
\item{CNOlist}{
The CNOlist used in optimization.
}
\item{indexList}{
An indexList, produced by indexFinder run on the model and the CNOlist above.
}
\item{plotPDF}{
TRUE or FALSE; for pdf output.
}
\item{tag}{
NULL or string; prefixes filenames with a tag (replaces
the default behavior).
}
\item{plotParams}{
A list of options related to the PDF and plotting outputs. (1) maxrow
is the maximum number of rows used to plot the results. See \link{plotOptimResultsPan}
for other fields.
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
\item{sizeFac}{
The scaling factor for the size term in the objective function, default to 0.0001.
}
\item{NAFac}{
The scaling factor for the NA term in the objective function, default to 1.
}
}
\value{
This function doesn't return anything, it only plots the graph in your graphic window.
}
\author{
A. MacNamara
}
\seealso{
plotOptimResultsPan
}
\examples{
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

# pre-process model
model <- preprocessing(CNOlistPB, modelPB)
indices = indexFinder(CNOlistPB, model)
fields4Sim <- prep4sim(model=model)
initBstring <- rep(1, length(model$reacID))

# optimize
opt1 <- gaBinaryDT(CNOlist=CNOlistPB, model=model, initBstring=initBstring,
verbose=TRUE, boolUpdates=10, maxTime=30, lowerB=0.8, upperB=10)

cutAndPlotResultsDT(
  model=model,
  CNOlist=CNOlistPB,
  bString=opt1$bString,
  plotPDF=FALSE,
  boolUpdates=10,
  lowerB=0.8, 
  upperB=10
)
}
