\name{cutAndPlot}
\alias{cutAndPlot}
\title{
Interface to cutAndPlotResultsT1 and T2 functions.
}
\description{
This function takes a model and cnolist as well as an optimised bitstring at T1
(and possibly T2). It then called cutAndPlotResultsT1 or T2.
}
\usage{
cutAndPlot(model, CNOlist, bStringT1, bStringT2, plotPDF =
FALSE, tag = NULL, show=TRUE,  tPt=2)
}
\arguments{
  \item{model}{
a model (the full one that was used for optimisation)
}
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{bStringT1}{
a bitstring for T1  as output by gaBinaryT1 (i.e. a vector of 1s and 0s)
}
  \item{bStringT1}{
a bitstring for T1  as output by gaBinaryT1 (i.e. a vector of 1s and 0s)
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
Time point index (default is 2). Note that T1 has index 2 and T2 has index 3.
}
}

\value{
This function doesn't return anything, it only plots the graph in your graphic window, and in a pdf file if asked
}
\author{
T. Cokelaer
}
\seealso{
cutAndPlotResultsT1
}
\examples{
#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model
results = preprocessing(CNOlistToy, ToyModel)
processedModel = results$model

cutAndPlot(processedModel,CNOlistToy,
	bStringT1=rep(1,length(processedModel$reacID)),	plotPDF=FALSE)
}
