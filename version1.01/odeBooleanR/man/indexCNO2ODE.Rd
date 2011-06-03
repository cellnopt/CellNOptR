\name{indexCNO2ODE}
\alias{indexCNO2ODE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{indexCNO2ODE
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
indexCNO2ODE(CNOlist, odePars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
}
  \item{odePars}{
%%     ~~Describe \code{odePars} here~~
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
function(CNOlist,odePars)
{ 
  CNOindexSignals=c();
  ODEindexSignals=c();
  for(i in 1:length(odePars$outputs))
  {
    tempIndex=which(CNOlist$namesSignals==odePars$outputs[i]);
    CNOindexSignals=c(CNOindexSignals,tempIndex);
    if(length(tempIndex)>0)ODEindexSignals=c(ODEindexSignals,i);
  }
  
  CNOindexStimuli=c();
  ODEindexStimuli=c();
  for(i in 1:length(odePars$inputs))
  {
    tempIndex=which(CNOlist$namesStimuli==odePars$inputs[i]);
    CNOindexStimuli=c(CNOindexStimuli,tempIndex);
    if(length(tempIndex)>0) ODEindexStimuli=c(ODEindexStimuli,i);
  }
  
  CNOindexInhibitors=c();
  ODEindexInhibitors=c();
  for(i in 1:length(odePars$outputs))
  {
    tempIndex=which(CNOlist$namesInhibitors==odePars$outputs[i]);
    CNOindexInhibitors=c(CNOindexInhibitors,tempIndex);
    if(length(tempIndex)>0)ODEindexInhibitors=c(ODEindexInhibitors,i);
  }
  res=list(CNOindexSignals=CNOindexSignals,ODEindexSignals=ODEindexSignals,
    CNOindexStimuli=CNOindexStimuli,ODEindexStimuli=ODEindexStimuli,
      CNOindexInhibitors=CNOindexInhibitors,ODEindexInhibitors=ODEindexInhibitors);
  
  return(res);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
