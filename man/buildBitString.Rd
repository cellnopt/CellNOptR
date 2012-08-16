\name{buildBitString}
\alias{buildBitString}
\title{
build final bit strings from a list of optimised bit strings at different
times.
}
\description{
This function takes as input a list of optimised bit strings at different time 
points (as returned by gaBinary functions). The input list must contain at least
one bit string that has the same length as the model to optimise. 

This function returns 2 lists of same length as the first bit string. The first
one is a map of all bit strings onto the first one. The second one is indicates
the time point at which each bit was optimised.

This function is used internally and should not be used by a user in principle.
However, it may be useful for post processing.

}
\usage{
buildBitString(bitStrings)
}
\arguments{
    \item{bitStrings}{
        a list of bit strings as returned by the optimisation at differnt time points.
    }
}
\value{
a list of 2 lists. The first one adds the input bitstrings together. The second
list (of same length as the first optimised bit string) indicates the time point
at which each bit was selected durint the optimisation processs.
}
\author{
T. Cokelaer
}
\examples{
#load data

res = buildBitString(c(1,1,0,1,0,0), c(0,1,0), c(1,0))
# res$bs = c(1,1,1,1,1,0)
# res$bsTimes = c(1,1,3,1,2,0)

}
