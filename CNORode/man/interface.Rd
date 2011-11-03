\name{interface}
\alias{interface}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
interface
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interface(cnolist, sif, indices, odeParameters, time = 1, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnolist}{
%%     ~~Describe \code{cnolist} here~~
}
  \item{sif}{
%%     ~~Describe \code{sif} here~~
}
  \item{indices}{
%%     ~~Describe \code{indices} here~~
}
  \item{odeParameters}{
%%     ~~Describe \code{odeParameters} here~~
}
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function(cnolist, sif, indices, odeParameters, time=1,verbose=0)
{	 
	interMat <- as.integer(as.vector(t(sif$interMat)))
	notMat <- as.integer(as.vector(t(sif$notMat)))
	nRows <- as.integer(dim(sif$interMat)[1])
	nCols <- as.integer(dim(sif$interMat)[2])
	verbose=as.integer(verbose);
	# ode 
	nPars <- as.integer(length(odeParameters))
	
	# cnolist
	timeSignals <- as.double(cnolist$timeSignals)
	valueInhibitors <- as.double(t(cnolist$valueInhibitors))
	valueSignals <- as.double(as.vector(cnolist$valueSignals[[time]]))
	valueStimuli <- as.double(t(cnolist$valueStimuli))
	nTimes = as.integer(length(cnolist$timeSignals))
	
	# [[1]] allows to access to the first object in the list and retrieve its dimensions
	nExperiments = as.integer(dim(cnolist$valueSignals[[time]]))
	
	#indices
	nSignals <- as.integer(length(indices$signals))
	
	indexSignals <- as.integer(as.vector(indices$signals))
	nStimuli <- as.integer(length(indices$stimulated))
	indexStimuli <- as.integer(as.vector(indices$stimulated))
	nInhibitors <- as.integer(length(indices$inhibited))
	indexInhibitors <- as.integer(as.vector(indices$inhibited))
	
	res = .Call("sim_logic_ode",interMat,notMat,nRows,nCols,nPars,timeSignals,
			valueInhibitors,valueSignals,valueStimuli,nTimes,nExperiments,nSignals,
			indexSignals,nStimuli,indexStimuli,nInhibitors,indexInhibitors,
			odeParameters,verbose);
	
	return(res);
	
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
