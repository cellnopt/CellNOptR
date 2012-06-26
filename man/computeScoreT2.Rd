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
computeScoreT2(CNOlist, model, simList, indexList, simResT1, bStringT1,
    bStringT2, sizeFac=0.0001, NAFac=1)
}


\arguments{
    \item{CNOlist}{
        a CNOlist structure, as created by \code{makeCNOlist}.
    }
    \item{model}{
        a model structure, as created by code{readSIF}, normally pre-processed 
        but that is not a requirement of this function.
    }
    \item{simList}{
        a simList as created by \code{prep4sim}, that has also already been cut to 
        contain only the reactions to be evaluated.
    }
    \item{indexList}{
        a list of indexes of the species stimulated/inhibited/measured in the 
        model, as created by \code{indexFinder}.
    }
    \item{simResT1}{
        Results of T1 step.
    }
    \item{bStringT1}{
        the best bitString at time T1
    }
    \item{bStringT2}{
        the bitString to use to compute the score at time T2
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
    data(CNOlistToy2,package="CellNOptR")
    data(ToyModel2,package="CellNOptR")
    res <- preprocessing(CNOlistToy2, ToyModel2)
    simlist <-prep4sim(res$model)
    bStringT1 = c(0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1)
    simT1<-simulateT1(CNOlist=CNOlistToy2, model=res$model,
            bStringT1=bStringT1, simList=simlist, indexList=res$indices)


    score = computeScoreT2(CNOlistToy2, res$model ,simlist, res$indices, simT1,
        bStringT1, bStringT2=c(1,0,1,0))
}



