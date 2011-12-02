\name{parEstimationLBodeSSm}
\alias{parEstimationLBodeSSm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parEstimationLBodeSSm(cnolist, model, ode_parameters = NULL, indices = NULL, maxeval = Inf, maxtime = 100, ndiverse = NULL, dim_refset = NULL, local_solver = NULL, time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50, nan_fac = 1)
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
  \item{maxeval}{
%%     ~~Describe \code{maxeval} here~~
}
  \item{maxtime}{
%%     ~~Describe \code{maxtime} here~~
}
  \item{ndiverse}{
%%     ~~Describe \code{ndiverse} here~~
}
  \item{dim_refset}{
%%     ~~Describe \code{dim_refset} here~~
}
  \item{local_solver}{
%%     ~~Describe \code{local_solver} here~~
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
  \item{nan_fac}{
%%     ~~Describe \code{nan_fac} here~~
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
    maxeval = Inf, maxtime = 100, ndiverse = NULL, dim_refset = NULL, 
    local_solver = NULL, time = 1, verbose = 0, transfer_function = 3, 
    reltol = 1e-04, atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, 
    maxErrTestsFails = 50, nan_fac = 1) 
{
    library(eSSmR)
    library(Rsolnp)
    if (is.null(ode_parameters)) {
        adjMat = incidence2Adjacency(model)
        ode_parameters = createLBodeContPars(adjMat, model$namesSpecies, 
            random = TRUE)
    }
    if (is.null(indices)) 
        indices <- indexFinder(cnolist, model, verbose = FALSE)
    problem = list()
    problem$f <- getLBodeContObjFunction(cnolist, model, ode_parameters, 
        indices, time, verbose, transfer_function, reltol, atol, 
        maxStepSize, maxNumSteps, maxErrTestsFails)
    problem$x_L <- ode_parameters$LB[ode_parameters$index_opt_pars]
    problem$x_U <- ode_parameters$UB[ode_parameters$index_opt_pars]
    problem$x_0 <- ode_parameters$parValues[ode_parameters$index_opt_pars]
    problem$int_var = 0
    problem$bin_var = 0
    opts = list()
    opts$maxeval = maxeval
    opts$maxtime = maxtime
    if (!is.null(local_solver)) 
        opts$local_solver = local_solver
    if (!is.null(ndiverse)) 
        opts$ndiverse = ndiverse
    if (!is.null(dim_refset)) 
        opts$dim_refset = dim_refset
    results = essR(problem, opts)
    ode_parameters$parValues[ode_parameters$index_opt_pars] = results$xbest
    ode_parameters$ssm_results = results
    return(ode_parameters)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
