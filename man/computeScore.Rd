\name{computeScore}
\alias{computeScore}

\title{
  Convert a model object in sif format
}

\description{
compute and return score of a model (cut using a bitstring).
}

\usage{
computeScore(CNOlist, Model, SimList, indexList, bitString)

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
a list that contains additional fields for the simulator, as created by prep4Sim
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


\details{
All links in the model are converted in sif format that is sourceNode-tab-sign-tab-targetNode. If there are ANDs, they are converted using dummy nodes called and# (e.g. A+B=C will be A-tab-1-tab-and1; B-tab-1-tab-and1; and1-tab-1-tab-C).
}

\value{
\item{score}{See gaBinaryT1 for details}
}

\author{
    A.Vatsiou, T. Cokelaer
}



