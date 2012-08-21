\name{cutAndPlotResultsTN}
\alias{cutAndPlotResultsTN}
\title{
Plot the results of an optimisation at tN
}
\description{
This function takes a model and an optimised bitstring, it cuts the model according to the bitstring and plots the results of the simulation along with the experimental data. This function is designed to work on results of a 2 step optimisation.
}


\usage{
cutAndPlotResultsTN(CNOlist, model, bStrings, plotPDF = FALSE, tag=NULL, show=TRUE)
}
\arguments{
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{model}{
a model (the full one that was used for optimisation)
}
  \item{bStrings}{
a list of bitstring at different time points 
}
  \item{plotPDF}{
TRUE or FALSE, tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
  \item{show}{
todo
}
}
\value{
This function doesn't return anything, it only plots the graph in your graphic window, and in a pdf file if asked
}
\author{
T. Cokelaer, A. MacNamara, Sarah Schrier, C. Terfve based on \link{cutAndPlotResultsT1}}

\seealso{
\link{gaBinaryT1}, \link{prep4sim}, \link{cutAndPlotResultsT1}
}
\examples{
#load data

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#pre-process model
results = preprocessing(CNOlistToy2, ToyModel2)
processedModel = results$model

#optimise t1
ToyT1opt2<-gaBinaryT1(
	CNOlist=CNOlistToy2,
	model=processedModel,
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

#Optimise T2
ToyT1opt2T2<-gaBinaryTN(
	CNOlist=CNOlistToy2,
	model=processedModel,
	bString=list(ToyT1opt2$bString),
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

cutAndPlotResultsTN(
	CNOlist=CNOlistToy2,
	model=processedModel,
	bStrings=list(ToyT1opt2$bString, bStringT2=ToyT1opt2T2$bString),
	plotPDF=FALSE)
}
