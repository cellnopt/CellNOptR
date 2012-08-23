\name{cutAndPlotResultsT1}
\alias{cutAndPlotResultsT1}
\title{
Plot the results of an optimisation at t1
}
\description{
This function takes a model and an optimised bitstring, it cuts the model according to the bitstring and plots the results of the simulation along with the experimental data.
}
\usage{
cutAndPlotResultsT1(model, bString, simList=NULL, CNOlist, indexList=NULL, plotPDF =
FALSE, tag = NULL, show=TRUE, tPt=CNOlist$timeSignals[2])
}

\arguments{
  \item{model}{
a model (the full one that was used for optimisation)
}
  \item{bString}{
a bitstring for T1  as output by gaBinaryT1 (i.e. a vector of 1s and 0s)
}
  \item{simList}{
deprecated argument kept for back compatibility. 
a simlist corresponding to the model, as output by prep4sim
}
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{indexList}{
deprecated argument kept for back compatibility. 
an indexList, produced by indexFinder ran on the model and the CNOlist above
}
  \item{plotPDF}{
TRUE or FALSE; tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
  \item{show}{
TRUE of FALSE; FALSE tells not to show the plot, which can be useful when scripting.
}
  \item{tPt}{
The number of time points in the data.}
}

\value{
This function returns a list of filenames that have been created.
}
\author{
C.Terfve, T. Cokelaer,A. MacNamara
}
\seealso{
\link{gaBinaryT1}
}
\examples{
#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model
model = preprocessing(CNOlistToy, ToyModel)

#optimise
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	model=model,
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

#plotting
cutAndPlotResultsT1(
	model=model,
	CNOlist=CNOlistToy,
	bString=ToyT1opt$bString,
	plotPDF=FALSE)
}
