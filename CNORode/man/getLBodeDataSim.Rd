\name{getLBodeDataSim}
\alias{getLBodeDataSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getLBodeDataSim(cnolist, model, ode_parameters = NULL, indices = NULL, time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnolist}{
%%     ~~Describe \code{cnolist} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{ode_parameters}{
%%     ~~Describe \code{ode_parameters} here~~
}
  \item{indices}{
%%     ~~Describe \code{indices} here~~
}
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{transfer_function}{
%%     ~~Describe \code{transfer_function} here~~
}
  \item{reltol}{
%%     ~~Describe \code{reltol} here~~
}
  \item{atol}{
%%     ~~Describe \code{atol} here~~
}
  \item{maxStepSize}{
%%     ~~Describe \code{maxStepSize} here~~
}
  \item{maxNumSteps}{
%%     ~~Describe \code{maxNumSteps} here~~
}
  \item{maxErrTestsFails}{
%%     ~~Describe \code{maxErrTestsFails} here~~
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
function (cnolist, model, ode_parameters = NULL, indices = NULL, 
    time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, 
    atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50) 
{
    adjMat = incidence2Adjacency(model)
    if (is.null(indices)) 
        indices <- indexFinder(cnolist, model, verbose = FALSE)
    if (is.null(ode_parameters)) 
        ode_parameters = createLBodeContPars(adjMat, model$namesSpecies)
    sim_function = getLBodeSimFunction(cnolist, model, adjMat, 
        indices1 = indices, odeParameters1 = ode_parameters$parValues, 
        time1 = time, verbose1 = verbose, transfer_function1 = transfer_function, 
        reltol1 = reltol, atol1 = atol, maxStepSize1 = maxStepSize, 
        maxNumSteps1 = maxNumSteps, maxErrTestsFails1 = maxErrTestsFails)
    sim = sim_function(cnolist, model, ode_parameters$parValues)
    sim = lapply(sim, function(x) x[, indices$signals])
    return(sim)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
