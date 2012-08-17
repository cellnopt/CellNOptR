\name{buildBitString}
\alias{buildBitString}
\title{
Concatenate final bit string and times of optimiation from a list of optimised bit strings at different
times.
}
\description{
This function takes as input a list of vectors (can be only one). Each vector
represent an optimised bit strings at different time points (as returned by 
the gaBinary functions). The first optimised bit string must represent your
model (i.e., have the same length). The length of the following vectors must
correspond to the number of zeros in the previous bitstring. For instance,
the following list of bit strings \code{bStrings = list(c(1,1,1,0,0,0),
c(1,0,0))} is correct whereas \code{bStrings = list(c(1,1,1,0,0,0), c(1,0))} is
incorrect.


This function is used internally \link{computeScoreTN} and
\link{simulateTN}. It should not be used by a user in principle.
However, it may be useful for post processing.

}
\usage{
buildBitString(bitStrings)
}
\arguments{
    \item{bitStrings}{
        a list of bit strings as returned by the optimisation at different time points.
    }
}


\value{
This function returns 2 vectors of same length as the first input vector. 
The first vector maps each vector on top of the first one. The second vector 
indicates the time point at which each bit was optimised (see example).
}
\author{
T. Cokelaer
}
\examples{

# Considering the optimised bitstrings at T1, T2 and T3 to be c(1,1,0,1,0,0),
# c(0,1,0) and c(0,1), then we can build the overall bitStrings as follows:

res = buildBitString(list(c(1,1,0,1,0,0), c(0,1,0), c(0,1)))

# The results bit string is accessed through the bs field:
res$bs
#[1] 1 1 0 1 1 1

# and times at which each bit is activated with the bsTimes field:
res$bsTimes

# res$bsTimes = c(1,1,0,1,2,3)

}
