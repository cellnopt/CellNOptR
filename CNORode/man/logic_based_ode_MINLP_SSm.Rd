\name{logic_based_ode_MINLP_SSm}
\alias{logic_based_ode_MINLP_SSm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logic_based_ode_MINLP_SSm(cnolist, model, ode_parameters = NULL, maxeval = Inf, maxtime = 100, ndiverse = NULL, dim_refset = NULL)
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
    maxtime = 100, ndiverse = NULL, dim_refset = NULL) 
{
    library(eSSmR)
    library(Rsolnp)
    if (is.null(ode_parameters)) {
        adjMat = incidence2Adjacency(model)
        ode_parameters = makeParameterList(adjMat, model$namesSpecies, 
            random = TRUE)
    }
    n_cont = length(ode_parameters$index_opt_pars)
    n_int = dim(model$interMat)[2]
    problem = list()
    problem$f <- get_logic_ode_MINLP_objective_function(cnolist, 
        model, ode_parameters, indices)
    problem$x_L <- c(ode_parameters$LB[ode_parameters$index_opt_pars], 
        matrix(0, 1, n_int))
    problem$x_U <- c(ode_parameters$UB[ode_parameters$index_opt_pars], 
        matrix(1, 1, n_int))
    problem$x_0 <- c(ode_parameters$parValues[ode_parameters$index_opt_pars], 
        as.integer(round(runif(12))))
    problem$int_var = 0
    problem$bin_var = n_int
    opts = list()
    opts$maxeval = maxeval
    opts$maxtime = maxtime
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
