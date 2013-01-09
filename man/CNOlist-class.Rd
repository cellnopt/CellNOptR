\name{CNOlist-class}
\docType{class}
\alias{CNOlist-class}
\alias{CNOlist,ANY}
\alias{CNOlist}
\alias{show,CNOlist-method}
\alias{length,CNOlist-method}
\alias{plot,CNOlist-method}
\alias{plot,CNOlist}
\alias{plot,CNOlist,CNOlist-method}

\alias{randomize,CNOlist-method}

\title{Class \code{"CNOlist"}}
\description{
This function takes as input the filename of a MIDAS file (or the list returned
by \link{makeCNOlist}) and returns an instance of CNOlist class. It provides an 
object oriented approach to manipulate CNOlist. This function calls
\link{readMIDAS} and \link{makeCNOlist}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNOlist", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{cues}:}{Object of class \code{"matrix"} }
    \item{\code{inhibitors}:}{Object of class \code{"matrix"}  }
    \item{\code{stimuli}:}{Object of class \code{"matrix"} }
    \item{\code{signals}:}{Object of class \code{"list"}  }
    \item{\code{timepoints}:}{Object of class \code{"vector"} timepoints contained in the signals matrix. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="CNOlist")}: ... }
%    \item{plot}{\code{signature(x="CNOlist", y="missing")}: ... }
    \item{plot}{\code{signature(x="CNOlist")}: Please see the help page for the
        \code{plot.CNOlist} method for more details.}
    \item{plot}{\code{signature(x="CNOlist", y="CNOlist")}: Please see the help
        page of \code{\link{plot.CNOlist}} method for more details.}
    \item{length}{\code{signature(x="CNOlist")}: returns length of the CNOlist
        that is the number of time points}
    \item{randomize}{\code{signature(x="CNOlist")}: randomize a cnolist using
        the \code{\link{randomizeCNOlist}} function from \code{CellNOptR}. See
        randomize.CNOlist help page for more details.}

	 }
}


\author{
T. Cokelaer
}
\seealso{
\code{\link{randomizeCNOlist}}, \code{\link{plotCNOlist}}, \code{\link{plotCNOlist2}}
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
new_cnolist = randomize(cnolist)
length(cnolist)
}
\keyword{classes}
