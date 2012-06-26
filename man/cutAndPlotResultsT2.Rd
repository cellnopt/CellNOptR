\name{cutAndPlotResultsT2}
\alias{cutAndPlotResultsT2}
\title{
Plot the results of an optimisation at t2
}
\description{
This function takes a model and an optimised bitstring, it cuts the model according to the bitstring and plots the results of the simulation along with the experimental data. This function is designed to work on results of a 2 step optimisation.
}
\usage{
cutAndPlotResultsT2(model, bStringT1, bStringT2, simList, CNOlist, indexList,
plotPDF = FALSE, tag=NULL, tPt=CNOlist$timeSignals[2:3])
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
  \item{simList}{
a simlist corresponding to the model, as output by prep4sim
}
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{indexList}{
an indexList, produced by indexFinder ran on the model and the CNOlist above
}
  \item{plotPDF}{
TRUE or FALSE, tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
  \item{tPt}{
todo.
}
}

\value{
This function doesn't return anything, it only plots the graph in your graphic window, and in a pdf file if asked
}
\author{
C. Terfve
}
\seealso{
gaBinaryT1, prep4sim, cutAndPlotResultsT1
}
\examples{
#load data

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#pre-process model


#optimise t1
results = preprocessing(CNOlistToy2, ToyModel2)
processedModel = results$model
indices = results$indices

ToyFields4Sim2<-prep4sim(processedModel)
initBstring2<-rep(1,length(processedModel$reacID))
ToyT1opt2<-gaBinaryT1(
	CNOlist=CNOlistToy2,
	model=processedModel,
	simList=ToyFields4Sim2,
	indexList=indices,
	initBstring=initBstring2,
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

#Optimise T2
SimToyT12<-simulateT1(
	CNOlist=CNOlistToy2,
	model=processedModel,
	bStringT1=ToyT1opt2$bString,
	simList=ToyFields4Sim2,
	indexList=indices)

ToyT1opt2T2<-gaBinaryT2(
	CNOlist=CNOlistToy2,
	model=processedModel,
	simList=ToyFields4Sim2,
	indexList=indices,
	bStringT1=ToyT1opt2$bString,
	simResT1=SimToyT12,
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

cutAndPlotResultsT2(
	model=processedModel,
	bStringT1=ToyT1opt2$bString,
	bStringT2=ToyT1opt2T2$bString,
	simList=ToyFields4Sim2,
	CNOlist=CNOlistToy2,
	indexList=indices,
	plotPDF=FALSE)
}
