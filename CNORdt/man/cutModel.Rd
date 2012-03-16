\name{cutModel}
\alias{cutModel}
\title{
Cuts model according to provided bit string}
\description{
This function cuts the expanded model according to a bit string (either optimized when used in CNORdt or can be user-provided as a vector of 0s and 1s the same length as the number of hyperedges in the expanded model.)
}
\usage{
cutModel(Model, SimList, bitString)
}
\arguments{
\item{Model}{
a Model that only contains the reactions to be evaluated.
}
\item{SimList}{
a list that contains additional fields for the simulator, as created by prep4Sim applied to the model above.
}
\item{bitString}{
a vector of 0s and 1s signifying which hyperedges to cut.
}
}
\value{
This function returns a list with elements:
\item{ModelCut}{the cut model.}
\item{SimListCut}{the cut SimList.}
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
data(ToyModelPB, package="CNORdt")

# pre-process model
indexOrig <- indexFinder(CNOlist=CNOlist, Model=Model, verbose=TRUE)
indexNONC <- findNONC(Model=Model, indexes=indexOrig, verbose=TRUE)
ModelCut <- cutNONC(Model=Model, NONCindexes=indexNONC)
indexNONCcut <- indexFinder(CNOlist=CNOlist, Model=ModelCut)
ModelCutCompress <- compressModel(Model=ModelCut, indexes=indexNONCcut)
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, Model=ModelCutCompress)
ModelCutCompressExpand <- expandGates(Model=ModelCutCompress)
fields4Sim <- prep4Sim(Model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

bitString = sample(c(0,1), length(ModelCutCompressExpand$reacID), replace=TRUE)
dataCut = cutModel(Model=ModelCutCompressExpand, SimList=fields4Sim, bitString=bitString)
}