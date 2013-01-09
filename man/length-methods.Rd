\name{length-method}
\docType{methods}
\alias{length.CNOlist}
\alias{length,CNOlist,ANY-method}



\title{Returns length of \code{"CNOlist"} object - methods}
\description{
 Returns the number of time points
}

\section{Usage}{
    \code{signature(x="CNOlist")}
}

\author{
T.Cokelaer
}
\seealso{
\code{\link{CNOlist}}
}
\examples{
showClass("CNOlist")

data(CNOlistToyMMB, package="CellNOptR")
cnolist = CNOlistToyMMB
length(cnolist)
}
\keyword{methods}
\keyword{CNOlist}
