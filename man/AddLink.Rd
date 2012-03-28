\name{AddLink}
\alias{AddLink}

\title{
  TODO
}

\description{
  TODO
}

\usage{
    AddLink(CueDown, SigUp, Model, Sign)
}


\arguments{
  \item{CueDown}{TODO}
  \item{SigUp}{TODO}
  \item{Model}{TODO}
  \item{Sign}{TODO}
}
\details{
TODO
}

\value{
\item{Model}{TODO}


}
\author{
    F.Eduati
}

\seealso{
    \link{readMIDAS}, 
}

\note{TODO}

\examples{
    \dontrun{
      library(CellNOptR)
      data(CNOlistToy,package="CellNOptR")
      data(ToyModel,package="CellNOptR")
      res<-AddLink(ToyModel, cnolist=CNOlistToy, compressed=c("TRAF6", "p38"))
    }

    library(CellNOptR)
    data(CNOlistToy,package="CellNOptR")
    data(ToyModel,package="CellNOptR")
    res<-AddLink(ToyModel, cnolist=CNOlistToy, compressed=c("TRAF6", "p38"))
    
}

