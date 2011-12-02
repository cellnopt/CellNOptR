\name{getLBodeSimFunction}
\alias{getLBodeSimFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getLBodeSimFunction(cnolist1, model1, adjMatrix1, indices1, odeParameters1, time1 = 1, verbose1 = 0, transfer_function1 = 3, reltol1 = 1e-04, atol1 = 0.001, maxStepSize1 = Inf, maxNumSteps1 = 1e+05, maxErrTestsFails1 = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnolist1}{
%%     ~~Describe \code{cnolist1} here~~
}
  \item{model1}{
%%     ~~Describe \code{model1} here~~
}
  \item{adjMatrix1}{
%%     ~~Describe \code{adjMatrix1} here~~
}
  \item{indices1}{
%%     ~~Describe \code{indices1} here~~
}
  \item{odeParameters1}{
%%     ~~Describe \code{odeParameters1} here~~
}
  \item{time1}{
%%     ~~Describe \code{time1} here~~
}
  \item{verbose1}{
%%     ~~Describe \code{verbose1} here~~
}
  \item{transfer_function1}{
%%     ~~Describe \code{transfer_function1} here~~
}
  \item{reltol1}{
%%     ~~Describe \code{reltol1} here~~
}
  \item{atol1}{
%%     ~~Describe \code{atol1} here~~
}
  \item{maxStepSize1}{
%%     ~~Describe \code{maxStepSize1} here~~
}
  \item{maxNumSteps1}{
%%     ~~Describe \code{maxNumSteps1} here~~
}
  \item{maxErrTestsFails1}{
%%     ~~Describe \code{maxErrTestsFails1} here~~
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
function (cnolist1, model1, adjMatrix1, indices1, odeParameters1, 
    time1 = 1, verbose1 = 0, transfer_function1 = 3, reltol1 = 1e-04, 
    atol1 = 0.001, maxStepSize1 = Inf, maxNumSteps1 = 1e+05, 
    maxErrTestsFails1 = 50) 
{
    simulate_logic_based_ode_model <- function(cnolist, sif, 
        odeParameters, indices = indices1, adjMatrix = adjMatrix1, 
        time = time1, verbose = verbose1, transfer_function = transfer_function1, 
        reltol = reltol1, atol = atol1, maxStepSize = maxStepSize1, 
        maxNumSteps = maxNumSteps1, maxErrTestsFails = maxErrTestsFails1) {
        interMat <- as.integer(as.vector(t(sif$interMat)))
        notMat <- as.integer(as.vector(t(sif$notMat)))
        adjMatrix <- as.integer(as.vector(t(adjMatrix)))
        nRows <- as.integer(dim(sif$interMat)[1])
        nCols <- as.integer(dim(sif$interMat)[2])
        verbose = as.integer(verbose)
        nPars <- as.integer(length(odeParameters))
        timeSignals <- as.double(cnolist$timeSignals)
        valueInhibitors <- as.double(t(cnolist$valueInhibitors))
        valueSignals <- as.double(t(cnolist$valueSignals[[time]]))
        valueStimuli <- as.double(t(cnolist$valueStimuli))
        nTimes = as.integer(length(cnolist$timeSignals))
        nExperiments = as.integer(dim(cnolist$valueSignals[[time]]))
        nSignals <- as.integer(length(indices$signals))
        indexSignals <- as.integer(as.vector(indices$signals))
        nStimuli <- as.integer(length(indices$stimulated))
        indexStimuli <- as.integer(as.vector(indices$stimulated))
        nInhibitors <- as.integer(length(indices$inhibited))
        indexInhibitors <- as.integer(as.vector(indices$inhibited))
        transfer_function = as.integer(transfer_function)
        reltol = as.double(reltol)
        atol = as.double(atol)
        maxStepSize = as.double(maxStepSize)
        maxNumSteps = as.integer(maxNumSteps)
        maxErrTestsFails = as.integer(maxErrTestsFails)
        break_at_1st_fail = as.integer(0)
        res = .Call("sim_logic_ode", interMat, notMat, adjMatrix, 
            nRows, nCols, nPars, timeSignals, valueInhibitors, 
            valueSignals, valueStimuli, nTimes, nExperiments, 
            nSignals, indexSignals, nStimuli, indexStimuli, nInhibitors, 
            indexInhibitors, odeParameters, verbose, transfer_function, 
            reltol, atol, maxStepSize, maxNumSteps, maxErrTestsFails, 
            break_at_1st_fail)
        return(res)
    }
    return(simulate_logic_based_ode_model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
