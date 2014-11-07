\name{prep4simFuzzy}
\alias{prep4simFuzzy}
\title{
Prepare a model for simulation
}
\description{
Adds to the model some fields that are used by the simulation engine and calls
prep4sim function from CellNOptR package
}
\usage{
prep4simFuzzy(model, paramsList, verbose=TRUE)
}
\arguments{
  \item{model}{
A model, as output by readSIF, normally pre-processed but that is not a
requirement of this function.
}
  \item{paramsList}{
A parameters list that must contain at least the CNOlist parameter in the
\emph{Data} field (param = list(data=cnolist)) and possibly optional arguments
related to the fuzzy logic (see \link{defaultParametersFuzzy})
}

 \item{verbose}{
a verbose option set to TRUE by default.
}


}
\details{
This adds fields that are necessary for the simulation engine for both Boolean and constrained Fuzzy logic simulation.
}
\value{
In addition to the fields returned by the prep4sim function of CellNOptR, this
function appends the following fields:
\item{finalCube}{stores, for each reac(row) the location of its inputs (col)}
\item{ixNeg}{stores, for each reac(row) and each input (col) whether it is a negative input}
\item{ignoreCube}{logical matrix of the same size as the 2 above, that tells whether the particular cell is filled or not}
\item{maxIx}{row vector that stores, for each reac, the location of its output}
\item{modelname}{stores the name of the model from which these fields were derived}
}
\author{
C. Terfve, modified by M.K. Morris and T. Cokelaer
}
\seealso{
\link{simFuzzyT1}
}
\examples{

    data(ToyModel, package="CellNOptR")
    data(CNOlistToy, package="CellNOptR")
    params <- defaultParametersFuzzy()
    params$data = CNOlistToy
    fields4sim <- prep4simFuzzy(ToyModel, params)
}
