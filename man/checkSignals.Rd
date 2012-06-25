\name{checkSignals}
\alias{checkSignals}
\title{
Check the CNOlist and model matching
}
\description{
This function checks that all the signals in the CNOlist match to species in the 
model. It also checks that the CNOlist and Model lists have the right format and 
contain the right fields. It is called by the \link{preprocessing} function so there is 
no need to call it directly anymore if you use the \link{preprocessing} function.
}
\usage{
checkSignals(CNOlist, Model)
}
\arguments{
    \item{CNOlist}{
        A CNOlist structure, as created by \link{makeCNOlist}.
    }
    \item{Model}{
        A model structure, as created by \link{readSIF}.
    }
}
\details{
If the formats are wrong, this function produces an error with an explanation
message. If the signals does not match the species, this function produces a 
warning that explains which signals don't match any species, and advises to 
remove them.
}

\author{
C. Terfve
}
\seealso{
\code{\link{makeCNOlist}}, \link{readSIF}, \link{preprocessing}
}
\examples{
data(CNOlistToy, package="CellNOptR")
data(ToyModel, package="CellNOptR")
checkSignals(CNOlistToy, ToyModel)
}


