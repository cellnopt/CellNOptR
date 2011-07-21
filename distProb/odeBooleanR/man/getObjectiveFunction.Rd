\name{getObjectiveFunction}
\alias{getObjectiveFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getObjectiveFunction(CNOlist1, odePars1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNOlist1}{
%%     ~~Describe \code{CNOlist1} here~~
}
  \item{odePars1}{
%%     ~~Describe \code{odePars1} here~~
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
function(CNOlist1,odePars1)
{
  evaluateModel<-function(x,CNOlist=CNOlist1,odePars=odePars1)
  { 
    odePars$paramValues[odePars$index_opt]=x;
    indexes=indexCNO2ODE(CNOlist,odePars); 
    simData=simulateData(CNOlist,indexes,odePars);
    times=CNOlist$timeSignals;
    expResults=CNOlist$valueSignals;  
    SimResults=convertTimeSeriesToList(simData,CNOlist);
    mySum=0;
    for(i in 1:length(SimResults))
    {
      mySum=mySum+(SimResults[[i]]-expResults[[i]])^2; 
    }
    mySum=sum(mySum);
    
    return(mySum)
  }
  return(evaluateModel);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
