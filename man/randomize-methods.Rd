\name{randomize-method}
\docType{methods}
\alias{randomize}
\alias{randomize.CNOlist}
\alias{randomize,CNOlist,ANY-method}



\title{plot a \code{"CNOlist"} object - methods}
\description{
 A randomize generator on a \code{CNOlist}. 
}

\section{Usage}{
    \code{signature(x="CNOlist", sd=0.1, minValue=0, maxValue=1, mode="gaussian")}: calls the
        \code{randomizeCNOlist} function.
}
\section{arguments}{
    \describe{
	\item{x}{The \code{CNOlist} object to plot}
	\item{sd}{The \code{CNOlist} object to plot}
	\item{minValue}{The \code{CNOlist} object to plot}
	\item{maxValue}{The \code{CNOlist} object to plot}
	\item{mode}{The \code{CNOlist} object to plot}
    }
}

\author{
T.Cokelaer
}
\seealso{
\code{\link{CNOlist}}, \code{\link{randomizeCNOlist}}
}
\examples{
showClass("CNOlist")

data(CNOlistToyMMB, package="CellNOptR")
cnolist = CNOlistToyMMB
# accessors:
# In version 1.3.30 and above, use the plot method instead of former plotCNOlist function.
new_cnolist = randomize(cnolist)
}
\keyword{methods}
\keyword{CNOlist}
