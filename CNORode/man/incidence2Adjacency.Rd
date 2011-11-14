\name{incidence2Adjacency}
\alias{incidence2Adjacency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
incidence2Adjacency
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
incidence2Adjacency(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model) 
{
    incidence = model$interMat
    nNodes = dim(incidence)[1]
    nEdges = dim(incidence)[2]
    adjacency = matrix(0, nNodes, nNodes)
    for (j in 1:nEdges) {
        for (i in 1:nNodes) {
            if (incidence[i, j] == 1) {
                node1 = i
                for (k in 1:nNodes) {
                  if (incidence[k, j] == -1) {
                    adjacency[k, node1] = 1
                  }
                }
            }
        }
    }
    return(adjacency)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
