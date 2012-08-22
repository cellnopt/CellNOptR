\name{cutAndPlotResultsT2}
\alias{cutAndPlotResultsT2}
\title{
Plot the results of an optimisation at t2 (deprecated since v1.3.28; use \link{cutAndPlotResultsTN} instead)}
\description{
This function takes a model and an optimised bitstring, it cuts the model according to the bitstring and plots the results of the simulation along with the experimental data. This function is designed to work on results of a 2 step optimisation.
}
\usage{
cutAndPlotResultsT2(model, bStringT1, bStringT2, CNOlist, simList=NULL, indexList=NULL,
plotPDF = FALSE, tag=NULL, show=TRUE, tPt=CNOlist$timeSignals[2:3])
}

\arguments{
  \item{model}{
a model (the full one that was used for optimisation)
}
  \item{bStringT1}{
a bitstring for T1 as output by gaBinaryT1 
}
  \item{bStringT2}{
a bitstring for T2 as output by gaBinaryT2
}
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{simList}{
deprecated argument kept for back compatibility.
a simlist corresponding to the model, as output by prep4sim
}
  \item{indexList}{
deprecated argument kept for back compatibility.
an indexList, produced by indexFinder ran on the model and the CNOlist above
}
  \item{plotPDF}{
TRUE or FALSE, tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
  \item{show}{
TRUE/FALSE 
}
  \item{tPt}{
todo.
}
}

\value{
This function doesn't return anything, it only plots the graph in your graphic window, and in a pdf file if asked
}
\author{
C. Terfve, T. Cokelaer, A. MacNamara}

\seealso{
\link{gaBinaryT1}, \link{cutAndPlotResultsT1}
}
\examples{
#load data

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#pre-process model
model = preprocessing(CNOlistToy2, ToyModel2)

#optimise t1
ToyT1<-gaBinaryT1(
	CNOlist=CNOlistToy2,
	model=model,
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

#Optimise T2

ToyT2<-gaBinaryT2(
	CNOlist=CNOlistToy2,
	model=model,
	bStringT1=ToyT1$bString,
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

cutAndPlotResultsT2(
	CNOlist=CNOlistToy2,
	model=model,
	bStringT1=ToyT1$bString,
	bStringT2=ToyT2$bString,
	plotPDF=FALSE)
}
