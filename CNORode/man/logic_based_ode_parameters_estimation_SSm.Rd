\name{logic_based_ode_parameters_estimation_SSm}
\alias{logic_based_ode_parameters_estimation_SSm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
logic_based_ode_parameters_estimation_SSm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logic_based_ode_parameters_estimation_SSm(cnolist, model, ode_parameters = NULL, maxeval = Inf, maxtime = 100, ndiverse = NULL, dim_refset = NULL, local_solver = NULL)
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
function (cnolist, model, ode_parameters = NULL, maxeval = Inf, 
    maxtime = 100, ndiverse = NULL, dim_refset = NULL, local_solver = NULL) 
{
    library(eSSmR)
    library(Rsolnp)
    if (is.null(ode_parameters)) {
        adjMat = incidence2Adjacency(model)
        ode_parameters = makeParameterList(adjMat, model$namesSpecies, 
            random = TRUE)
    }
    problem = list()
    problem$f <- get_logic_ode_continuous_objective_function(cnolist, 
        model, ode_parameters, indices)
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
    return(essR(problem, opts))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
