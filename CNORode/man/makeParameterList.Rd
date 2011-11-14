\name{makeParameterList}
\alias{makeParameterList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
makeParameterList
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeParameterList(adjMat, namesSpecies, LB_n = 1, LB_k = 0.1, LB_tau = 0.01, UB_n = 5, UB_k = 0.9, UB_tau = 10, default_n = 3, default_k = 0.5, default_tau = 1, LB_in = c(), UB_in = c(), opt_n = TRUE, opt_k = TRUE, opt_tau = TRUE, random = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjMat}{
%%     ~~Describe \code{adjMat} here~~
}
  \item{namesSpecies}{
%%     ~~Describe \code{namesSpecies} here~~
}
  \item{LB_n}{
%%     ~~Describe \code{LB_n} here~~
}
  \item{LB_k}{
%%     ~~Describe \code{LB_k} here~~
}
  \item{LB_tau}{
%%     ~~Describe \code{LB_tau} here~~
}
  \item{UB_n}{
%%     ~~Describe \code{UB_n} here~~
}
  \item{UB_k}{
%%     ~~Describe \code{UB_k} here~~
}
  \item{UB_tau}{
%%     ~~Describe \code{UB_tau} here~~
}
  \item{default_n}{
%%     ~~Describe \code{default_n} here~~
}
  \item{default_k}{
%%     ~~Describe \code{default_k} here~~
}
  \item{default_tau}{
%%     ~~Describe \code{default_tau} here~~
}
  \item{LB_in}{
%%     ~~Describe \code{LB_in} here~~
}
  \item{UB_in}{
%%     ~~Describe \code{UB_in} here~~
}
  \item{opt_n}{
%%     ~~Describe \code{opt_n} here~~
}
  \item{opt_k}{
%%     ~~Describe \code{opt_k} here~~
}
  \item{opt_tau}{
%%     ~~Describe \code{opt_tau} here~~
}
  \item{random}{
%%     ~~Describe \code{random} here~~
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
function (adjMat, namesSpecies, LB_n = 1, LB_k = 0.1, LB_tau = 0.01, 
    UB_n = 5, UB_k = 0.9, UB_tau = 10, default_n = 3, default_k = 0.5, 
    default_tau = 1, LB_in = c(), UB_in = c(), opt_n = TRUE, 
    opt_k = TRUE, opt_tau = TRUE, random = FALSE) 
{
    isState = getStates(adjMat)
    numStates = as.integer(sqrt(length(adjMat)))
    parNames = c()
    parValues = c()
    index_k = c()
    index_n = c()
    index_tau = c()
    LB = c()
    UB = c()
    index_opt_pars = c()
    count = 0
    for (j in 1:numStates) {
        if (isState[j]) {
            jCol = adjMat[, j]
            inputs = which(as.logical(jCol))
            numInputs = length(inputs)
            for (i in 1:numInputs) {
                count = count + 1
                parNames[count] = paste(namesSpecies[inputs[i]], 
                  "_n_", namesSpecies[j], sep = "")
                parValues[count] = default_n
                index_k = c(index_k, count)
                LB[count] = LB_n
                UB[count] = UB_n
                if (opt_n) 
                  index_opt_pars = c(index_opt_pars, count)
                count = count + 1
                parNames[count] = paste(namesSpecies[inputs[i]], 
                  "_k_", namesSpecies[j], sep = "")
                parValues[count] = default_k
                index_n = c(index_k, count)
                LB[count] = LB_k
                UB[count] = UB_k
                if (opt_k) 
                  index_opt_pars = c(index_opt_pars, count)
            }
            count = count + 1
            parNames[count] = paste("tau_", namesSpecies[j], 
                sep = "")
            parValues[count] = default_tau
            index_tau = c(index_tau, count)
            LB[count] = LB_tau
            UB[count] = UB_tau
            if (opt_tau) 
                index_opt_pars = c(index_opt_pars, count)
        }
    }
    if (length(LB_in) == length(LB)) 
        LB = LB_in
    if (length(UB_in) == length(UB)) 
        UB = UB_in
    if (random) 
        parValues = LB[index_opt_pars] + ((UB[index_opt_pars] - 
            LB[index_opt_pars]) * runif(length(index_opt_pars)))
    parList = list(parNames = parNames, parValues = parValues, 
        index_opt_pars = index_opt_pars, index_n = index_n, index_k = index_k, 
        index_tau = index_tau, LB = LB, UB = UB)
    return(parList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
