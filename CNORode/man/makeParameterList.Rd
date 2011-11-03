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
makeParameterList(adjMat, namesSpecies)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjMat}{
%%     ~~Describe \code{adjMat} here~~
}
  \item{namesSpecies}{
%%     ~~Describe \code{namesSpecies} here~~
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
function(adjMat,namesSpecies)
{
	isState=getStates(adjMat);
	numStates=as.integer(sqrt(length(adjMat)));
	parNames=c();
	parValues=c();
	count=0;
	for(j in 1:numStates)
	{ 
		if(isState[j])
		{
			jCol=adjMat[,j];
			inputs=which(as.logical(jCol));
			numInputs=length(inputs);
			for(i in 1:numInputs)
			{
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_n_",namesSpecies[j],sep="");
				parValues[count]=1;
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_k_",namesSpecies[j],sep="");
				parValues[count]=0.5;
			}
			count=count+1;
			parNames[count]=paste("tau_",namesSpecies[j],sep="");
			parValues[count]=1;
		}
	}
	parList=list(parNames=parNames,parValues=parValues);
	
	return(parList);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
