\name{findNONC}
\alias{findNONC}
\title{
Find the indexes of the non-observable and non controllable species
}
\description{
This function finds the indexes of the non-observable and non controllable species and returns the indices, in the model, of the species to remove
}
\usage{
findNONC(Model, indexes, verbose=FALSE)
}
\arguments{
    \item{Model}{
        a model structure, as created by \link{readSIF}
    }
    \item{indexes}{
        a list of indexes of the species stimulated/inhibited/measured, as created by \code{indexFinder} from a model.        
    }
    \item{verbose}{
        verbose option (default to FALSE)
    }
}
\details{
This function uses the function floyd.warshall.all.pairs.sp from the package RBGL.  Non observable nodes are those that do not have a path to any measured species in the model, whereas non controllable nodes are those that do not receive any information from a species that is perturbed in the data.
}
\value{
a vector of indices of species to remove
}
\author{
C. Terfve
}
\seealso{
cutNONC, indexFinder, readSIF
}
\examples{
    data(CNOlistToy, package="CellNOptR")
    data(ToyModel, package="CellNOptR")
    checkSignals(CNOlistToy, ToyModel)
    indicesToy <- indexFinder(CNOlistToy, ToyModel)
    ToyNCNOindices <- findNONC(ToyModel, indicesToy)
}
