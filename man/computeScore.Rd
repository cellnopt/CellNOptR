\name{computeScore}
\alias{computeScore}

\title{
  Compute the score of a model/data set using a bitString to cut the model.
}

\description{
See \code{gaBinaryT1}
}

\usage{
computeScore(CNOlist, Model, SimList, indexList, bitString, sizeFac=0.0001,
NAFac=1)

}


\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at
t1)
}
  \item{Model}{
a Model list
}
  \item{SimList}{
a list that contains additional fields for the simulator, as created by prep4sim
applied to the model above
}
  \item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by
indexfinder applied on the model and CNOlist above
}
  \item{bitString}{
a bistring to select a sub model
}
}



\value{
\item{score}{See gaBinaryT1 for details}
}

\author{
    T. Cokelaer
}



