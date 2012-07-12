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
a Model that only contains the reactions to be evaluated.
}
\item{simList}{
a list that contains additional fields for the simulator, as created by prep4sim applied to the model above.
}
\item{bitString}{
a vector of 0s and 1s signifying which hyperedges to cut.
}
}
\value{
This function returns a list with elements:
\item{modelCut}{the cut model.}
\item{simListCut}{the cut SimList.}
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
indexOrig <- indexFinder(CNOlist=CNOlist, model=Model, verbose=TRUE)
indexNONC <- findNONC(model=Model, indexes=indexOrig, verbose=TRUE)
ModelCut <- cutNONC(model=Model, NONCindexes=indexNONC)
indexNONCcut <- indexFinder(CNOlist=CNOlist, model=ModelCut)
ModelCutCompress <- compressModel(model=ModelCut, indexes=indexNONCcut)
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, model=ModelCutCompress)
ModelCutCompressExpand <- expandGates(model=ModelCutCompress)
fields4Sim <- prep4sim(model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

bitString = sample(c(0,1), length(ModelCutCompressExpand$reacID), replace=TRUE)
dataCut = cutModel(model=ModelCutCompressExpand, simList=fields4Sim, bitString=bitString)
}
