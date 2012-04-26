\name{sif2graph}
\alias{sif2graph}

\title{
  Convert sif to graph
}

\description{
  This function converts a network form sif format to graph format.
}

\usage{
sif2graph(sif)
}


\arguments{
  \item{sif}{the name of a sif file}
}
\details{
This function takes in a single argument, sifFile, that points to a previous knowledge network in .sif format i.e. sourceNode-tab-sign-tab-targetNode. If there are ANDs they should be introduced as dummy nodes called and# (don't forget the number after "and" otherwise this won't be recognised).  Please be aware that "and" nodes are not expected to be negated, i.e. there are not supposed to be !and1=xyz because that amounts to inverting the sign of all inputs of and1, which is more simply done at the inputs level.
}

\value{
\item{g}{the corresponding graph}
}

\author{
    F.Eduati
}

\seealso{
\link{model2sif}
}


