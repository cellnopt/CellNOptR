\name{cutModel}
\alias{cutModel}
\title{
Cuts model according to provided bit string}
\description{
This function cuts the expanded model according to a bit string (either optimized when used in CNORdt or can be user-provided as a vector of 0s and 1s the same length as the number of hyperedges in the expanded model.)
}
\usage{
cutModel(model, simList, bitString)
}
\arguments{
\item{model}{
A Model that only contains the reactions to be evaluated.
}
\item{simList}{
A list that contains additional fields for the simulator, as created by prep4sim applied to the model above.
}
\item{bitString}{
A vector of 0s and 1s signifying which hyperedges to cut.
}
}
\value{
This function returns a list with elements:
\item{modelCut}{The cut model.}
\item{simListCut}{The cut SimList.}
}

\author{
A. MacNamara
}
\seealso{
gaBinaryTimeScale, getFitTimeScale
}
\examples{
# this function is normally not visible to the user. It can be used though with a user-provided bit string:

library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

# pre-process model
model <- preprocessing(CNOlistPB, modelPB)
fields4Sim <- prep4sim(model=model)
initBstring <- rep(1, length(model$reacID))

bitString = sample(c(0,1), length(model$reacID), replace=TRUE)
dataCut = cutModel(model=model, simList=fields4Sim, bitString=bitString)
}
