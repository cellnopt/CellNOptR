\name{computeScoreFuzzy}
\alias{computeScoreFuzzy}

\title{
  Compute Score of a model compared to the data for a given intString.
}

\description{
compute and return score of a model (cut using a bitstring).
}

\usage{
computeScoreFuzzy(CNOlist, model,simList=NULL, indexList=NULL, paramsList, 
    intString=NULL, sizeFac=0.0001,NAFac=1)
}


\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at
t1)
}
  \item{model}{
a model list
}
  \item{simList}{
a list that contains additional fields for the simulator, as created by prep4sim
applied to the model above
}
  \item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by
indexfinder applied on the model and CNOlist above
}
    \item{paramsList}{list of parameters. See \link{defaultParametersFuzzy}}
    \item{intString}{
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
\item{score}{See gaBinaryT1 for details}
}

\author{
    T. Cokelaer
}



