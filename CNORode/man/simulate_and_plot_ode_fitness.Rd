\name{simulate_and_plot_ode_fitness}
\alias{simulate_and_plot_ode_fitness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulate_and_plot_ode_fitness(cnolist, model, ode_parameters = NULL, indices = NULL, adjMatrix = NULL, time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50, plot_index_signals = NULL, plot_index_experiments = NULL, plot_index_cues = NULL)
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
  \item{adjMatrix}{
%%     ~~Describe \code{adjMatrix} here~~
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
  \item{plot_index_signals}{
%%     ~~Describe \code{plot_index_signals} here~~
}
  \item{plot_index_experiments}{
%%     ~~Describe \code{plot_index_experiments} here~~
}
  \item{plot_index_cues}{
%%     ~~Describe \code{plot_index_cues} here~~
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
function
(
		cnolist,				    model,					    ode_parameters=NULL,
		indices=NULL,			    adjMatrix=NULL,			 	time=1,
		verbose=0, 				    transfer_function=3,		reltol=1e-4,
		atol=1e-3,				    maxStepSize=Inf,		 	maxNumSteps=100000,
		maxErrTestsFails=50,		plot_index_signals=NULL,	plot_index_experiments=NULL,
		plot_index_cues=NULL
)
{
	if(is.null(plot_index_experiments))plot_index_experiments=1:dim(cnolist$valueCues)[1];
	if(is.null(plot_index_cues))plot_index_cues=1:dim(cnolist$valueCues)[2];
	if(is.null(plot_index_signals))plot_index_signals=1:dim(cnolist$valueSignals[[1]])[1];
	if(is.null(indices))indices=indexFinder(cnolist,model);
	if(is.null(adjMatrix))adjMatrix=incidence2Adjacency(model);
	if(is.null(ode_parameters))ode_parameters=makeParameterList(adjMatrix,model$namesSpecies);

	sim_data=get_logic_based_ode_data_simulation(cnolist,model,
			ode_parameters,indices,time,verbose,transfer_function,
			reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);

	times=cnolist$timeSignals;
	sim_data=lapply(sim_data,function(x)x[plot_index_experiments,plot_index_signals]);
	expResults=lapply(cnolist$valueSignals,function(x)x[plot_index_experiments,plot_index_signals]);

	namesSignals=cnolist$namesSignals[plot_index_signals];
	namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);

	valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
	valueCues=as.matrix(valueCues);
	valueCues[which(valueCues>0)]=1;
	valueCues=valueCues[plot_index_experiments,plot_index_cues];
	names(valueCues)=namesCues[plot_index_cues];

	plotOptimResults(SimResults=sim_data,expResults=expResults,
			times=times,namesCues=namesCues,namesSignals=namesSignals,
			valueCues=valueCues);
			
  return(sim_data);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
