\name{simulateModel}
\alias{simulateModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simulateModel
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulateModel(CNOlist, indexes, odePars, times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
}
  \item{indexes}{
%%     ~~Describe \code{indexes} here~~
}
  \item{odePars}{
%%     ~~Describe \code{odePars} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
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
function(CNOlist,indexes,odePars,times)
{
  numExperiments=length(CNOlist$valueStimuli[,1]);
  res=list();
  for(i in 1:numExperiments)
  {
    odePars$paramValues[odePars$index_inputs[indexes$ODEindexStimuli]]=
      CNOlist$valueStimuli[i,indexes$CNOindexStimuli];
    odePars$paramValues[odePars$index_inh[indexes$ODEindexInhibitors]]=
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
      y0=matrix(0,1,length(odePars$outputs));
    y0[indexes$ODEindexSignals]=CNOlist$valueSignals[[1]][i,indexes$CNOindexSignals]
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
    tempRes=createInSilicoData(odePars,y0,times,"myODE");
    res[[i]]=tempRes;    
  }
  
  return(res);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
