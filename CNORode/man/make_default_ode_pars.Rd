\name{make_default_ode_pars}
\alias{make_default_ode_pars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make_default_ode_pars
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_default_ode_pars(interMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interMat}{
%%     ~~Describe \code{interMat} here~~
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
function(interMat)
{
	adjMat=getAdjMat(interMat);
	n_nodes=dim(adjMat)[1];
	ode_parameters=c();
	
	for (j in 1:n_nodes) 
	{
		isState=FALSE;
		for (i in 1:n_nodes) 
		{
			if(adjMat[i,j])
			{
				ode_parameters=c(ode_parameters,3);
				ode_parameters=c(ode_parameters,0.5);
				isState=TRUE;
			}	
		}
		if(isState)ode_parameters=c(ode_parameters,1);
	}
	return(ode_parameters)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
