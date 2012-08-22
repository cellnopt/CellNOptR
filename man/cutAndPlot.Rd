\name{cutAndPlot}
\alias{cutAndPlot}
\title{
Interface to cutAndPlotResults functions.
}
\description{
This function takes a model and cnolist as well as a list of optimised bitstring at different time points. It calls the appropriate cutAndPlotResultsTX function.
}
\usage{
cutAndPlot(CNOlist, model,  bStrings, plotPDF=FALSE, tag=NULL, show=TRUE  )
}
\arguments{
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{model}{
a model (the full one that was used for optimisation)
}
  \item{bStrings}{
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
}

\value{
This function returns nothing. It plots a graph in your graphic window and sqve it in a file if asked
}
\author{
T. Cokelaer
}
\seealso{
\link{cutAndPlotResultsT1}
}
\examples{
#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model
model = preprocessing(CNOlistToy, ToyModel)

cutAndPlot(CNOlistToy, model,
    bStrings=list(rep(1,length(model$reacID))),	plotPDF=FALSE)
}
