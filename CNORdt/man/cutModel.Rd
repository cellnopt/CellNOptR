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
\item{CNOlist}{
a CNOlist.
}
\item{Model}{
a Model that only contains the reactions to be evaluated.
}
\item{bitString}{
a vector of 0s and 1s signifying which hyperedges to cut.
}
}
\details{
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
library(CNORdiscreteTime)
data(CNOlistToyPB, package="CNORdiscreteTime")
data(ToyModelPB, package="CNORdiscreteTime")

# pre-process model
indexOrig <- indexFinder(CNOlist=CNOlist, Model=Model, verbose=T)
indexNONC <- findNONC(Model=Model, indexes=indexOrig, verbose=T)
ModelCut <- cutNONC(Model=Model, NONCindexes=indexNONC)
indexNONCcut <- indexFinder(CNOlist=CNOlist, Model=ModelCut)
ModelCutCompress <- compressModel(Model=ModelCut, indexes=indexNONCcut)
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, Model=ModelCutCompress)
ModelCutCompressExpand <- expandGates(Model=ModelCutCompress)
fields4Sim <- prep4Sim(Model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

bitString = sample(c(0,1), length(ModelCutCompressExpand$reacID), replace=T)
dataCut = cutModel(Model=ModelCutCompressExpand, SimList=fields4Sim, bitString=bitString)
}