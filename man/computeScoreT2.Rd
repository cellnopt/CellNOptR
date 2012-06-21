\name{computeScoreT2}
\alias{computeScoreT2}

\title{
  Compute the score at T2 of a model/data set using a bitString to cut the model.
}

\description{
The bitString made of 0 and 1 allows to select a submodel from the model provided. 
Then, the simulator function are called to compute the objective function. The 
sizeFac and NAFac are penalties added to the final score as described in 
\code{gaBinaryT2}.
}

\usage{
computeScoreT2(CNOlist, Model, SimList, indexList, SimResT1, bitString, bStringTimes, sizeFac=0.0001, NAFac=1)
}


\arguments{
    \item{CNOlist}{
        a CNOlist structure, as created by \code{makeCNOlist}.
    }
    \item{Model}{
        a model structure, as created by code{readSIF}, normally pre-processed 
        but that is not a requirement of this function.
    }
    \item{SimList}{
        a SimList as created by \code{prep4sim}, that has also already been cut to 
        contain only the reactions to be evaluated.
    }
    \item{indexList}{
        a list of indexes of the species stimulated/inhibited/measured in the 
        model, as created by \code{indexFinder}.
    }
    \item{SimResT1}{
    Results of T1 step.
    }
    \item{bitString}{
        a bitstring of the same size as the number of reactions in the model above
    }
    \item{sizeFac}{
    the scaling factor for the size term in the objective function, default to 0.0001
    }
    \item{NAFac}{
    the scaling factor for the NA term in the objective function, default to 1    
    }
}



\value{
    \item{score}{See \code{gaBinaryT2} for details}
}

\author{
    T. Cokelaer
}

\examples{
    data(CNOlistToy,package="CellNOptR")
    data(ToyModel,package="CellNOptR")
    res <- preprocessing(CNOlistToy,ToyModel)
    simlist <-prep4sim(res$model)
    score = computeScoreT2(CNOlistToy, res$model,simlist, res$indices, c())
}



