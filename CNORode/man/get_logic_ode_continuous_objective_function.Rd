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
get_logic_ode_continuous_objective_function(cnolist1, model1, ode_parameters1, indices1, time1 = 1, verbose1 = 0, transfer_function1 = 3, reltol1 = 1e-04, atol1 = 0.001, maxStepSize1 = Inf, maxNumSteps1 = 1e+05, maxErrTestsFails1 = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnolist1}{
%%     ~~Describe \code{cnolist1} here~~
}
  \item{model1}{
%%     ~~Describe \code{model1} here~~
}
  \item{ode_parameters1}{
%%     ~~Describe \code{ode_parameters1} here~~
}
  \item{indices1}{
%%     ~~Describe \code{indices1} here~~
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
function(cnolist1,model1,ode_parameters1,indices1,
		time1=1,verbose1=0, transfer_function1=3,reltol1=1e-4,atol1=1e-3,maxStepSize1=Inf,maxNumSteps1=100000,
		maxErrTestsFails1=50)
{
	adjMatrix1=incidence2Adjacency(model1);
	sim_function1=get_simulation_function(cnolist1,model1,adjMatrix=adjMatrix1,
			indices=indices1, odeParameters=ode_parameters1$parValues, time=time1,verbose=verbose1,
			transfer_function=transfer_function1,reltol=reltol1,atol=atol1,maxStepSize=maxStepSize1,
			maxNumSteps=maxNumSteps1,maxErrTestsFails=maxErrTestsFails1)

	logic_ode_continuous_objective_function<-
			function(x,cnolist=cnolist1,model=model1,adjMatrix=adjMatrix1,
					indices=indices1,ode_parameters=ode_parameters1,sim_function=sim_function1)
	{
		ode_parameters$parValues[ode_parameters$index_opt_pars]=x;
		sim=sim_function(cnolist,model,indices,ode_parameters$parValues);
		sim<-unlist(lapply(sim,function(x) x[,indices$signals]));
		measured_values=unlist(cnolist$valueSignals);                                                    
		NaNs=which(is.na(sim));
		not_NaNs=which(!is.na(sim));
		error=sum((sim[not_NaNs]-measured_values[not_NaNs])^2);
		if(is.na(error))error=0;
		res=(error+length(NaNs))/length(sim);
		return(res);
	}
	return(logic_ode_continuous_objective_function);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
