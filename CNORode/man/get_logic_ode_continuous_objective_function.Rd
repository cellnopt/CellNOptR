\name{get_logic_ode_continuous_objective_function}
\alias{get_logic_ode_continuous_objective_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get_logic_ode_continuous_objective_function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_logic_ode_continuous_objective_function(cnolist, model, ode_parameters, indices, time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50)
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
function (cnolist, model, ode_parameters, indices, time = 1, 
    verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, 
    maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50) 
{
    adjMatrix = incidence2Adjacency(model)
    sim_function = get_simulation_function(cnolist, model, adjMatrix1 = adjMatrix, 
        indices1 = indices, odeParameters1 = ode_parameters$parValues, 
        time1 = time, verbose1 = verbose, transfer_function1 = transfer_function, 
        reltol1 = reltol, atol1 = atol, maxStepSize1 = maxStepSize, 
        maxNumSteps1 = maxNumSteps, maxErrTestsFails1 = maxErrTestsFails)
    logic_ode_continuous_objective_function <- function(x, cnolist1 = cnolist, 
        model1 = model, adjMatrix1 = adjMatrix, indices1 = indices, 
        ode_parameters1 = ode_parameters, sim_function1 = sim_function) {
        ode_parameters1$parValues[ode_parameters1$index_opt_pars] = x
        sim = sim_function1(cnolist1, model1, ode_parameters1$parValues)
        sim <- unlist(lapply(sim, function(x) x[, indices1$signals]))
        measured_values = unlist(cnolist1$valueSignals)
        NaNs = which(is.na(sim))
        not_NaNs = which(!is.na(sim))
        error = sum((sim[not_NaNs] - measured_values[not_NaNs])^2)
        if (is.na(error)) 
            error = 0
        res = (error + length(NaNs))/length(sim)
        return(res)
    }
    return(logic_ode_continuous_objective_function)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
