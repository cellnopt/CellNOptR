\name{plot-method}
\docType{methods}
\alias{plot.CNOlist}
\alias{plot,CNOlist,ANY-method}
\alias{plot,CNOlist,ANY,ANY-method}



\title{plot a \code{"CNOlist"} object - methods}
\description{
 A plot method for \code{CNOlist}. 
}

\section{Usage}{
    \code{signature(x="CNOlist")}: If only one argument (a CNOlist) is provided,
        this method calls \code{\link{plotCNOlist}}.

    \code{signature(x="CNOlist", y="CNOlist")}: If 2 arguments (2 CNOlist
        objects) are provided, this method calls \code{\link{plotCNOlist2}}.
}
\section{arguments}{
    \describe{
	\item{x}{The \code{CNOlist} object to plot}
	\item{y}{an optional \code{CNOlist} object to plot}
    }
}


\author{
T. Cokelaer
}
\seealso{
\code{\link{CNOlist}}
}
\examples{
showClass("CNOlist")

files<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
cnolist  = CNOlist(files[[1]])
# accessors:
cnolist@cues
cnolist@inhibitors
cnolist@signals
cnolist@timepoints
cnolist@stimuli
# In version 1.3.30 and above, use the plot method instead of former plotCNOlist function.
plot(cnolist)
}
\keyword{methods}
\keyword{CNOlist}
