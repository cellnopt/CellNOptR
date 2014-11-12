\name{computeScoreDT}
\alias{computeScoreDT}
\title{
Compute the score of a model/data (time-course) set using a bitString to cut the model.
}

\description{
The bitString made of 0s and 1s encodes a submodel from the model provided. 
Then, the simulator function is called to compute the objective function. The 
sizeFac and NAFac are penalties added to the final score as described in 
\code{gaBinaryDT}. The indexList and simList arguments can be provided to speed up the code, otherwise 
they are recomputed from the CNOlist and model.
}

\usage{
computeScoreDT(CNOlist, model, bString, simList=NULL, indexList=NULL,
sizeFac=0.0001, NAFac=1, boolUpdates, lowerB=lowerB, upperB=upperB)
}


\arguments{
\item{CNOlist}{
A CNOlist structure, as created by \code{makeCNOlist}.
}
\item{model}{
A model structure, as created by code{readSIF}, normally pre-processed 
but that is not a requirement of this function.
}
\item{bString}{
A bit string of the same size as the number of reactions in the model above.
}
\item{simList}{
If provided, simList should be created by \code{prep4sim}, and has also already been cut to 
contain only the reactions to be evaluated. 
}
\item{indexList}{
If provided, indexList should contain a list of indexes of the species stimulated/inhibited/measured in the 
model, as created by \code{indexFinder}.
}
\item{sizeFac}{
The scaling factor for the size term in the objective function, default to 0.0001.
}
\item{NAFac}{
The scaling factor for the NA term in the objective function, default to 1.
}
\item{boolUpdates}{
The number of synchronous updates performed by the boolean simulator.
}
\item{lowerB}{
The lower bound for the optimized value of the scaling factor.
}
\item{upperB}{
The upper bound for the optimized value of the scaling factor.
}
}

\value{
\item{score}{See \code{gaBinaryT1} for details.}
}

\author{
A. MacNamara
}

\examples{
library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

# pre-process model
model = preprocessing(CNOlistPB, modelPB)

# compute score
score = computeScoreDT(CNOlistPB, model, bString=rep(1,16),
boolUpdates=10, lowerB=0.8, upperB=10)
}
