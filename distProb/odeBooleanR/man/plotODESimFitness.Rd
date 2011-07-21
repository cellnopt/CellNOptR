\name{plotODESimFitness}
\alias{plotODESimFitness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotODESimFitness(SimResults, CNOList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SimResults}{
%%     ~~Describe \code{SimResults} here~~
}
  \item{CNOList}{
%%     ~~Describe \code{CNOList} here~~
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
function(SimResults,CNOList)
{ 
  times=CNOlist$timeSignals;
  expResults=CNOlist$valueSignals;  
  SimResults=convertTimeSeriesToList(SimResults,CNOlist);
  
  namesSignals=CNOlist$namesSignals;
  
  namesCues=c(CNOlist$namesStimuli,CNOlist$namesInhibitors);
  
  valueCues=cbind(CNOlist$valueStimuli,CNOlist$valueInhibitors)
  valueCues=as.matrix(valueCues);
  valueCues[which(valueCues>0)]=1;
  names(valueCues)=namesCues;
  
  plotOptimResults(SimResults=SimResults,expResults=expResults,times=times,namesCues=namesCues,namesSignals=namesSignals,valueCues=valueCues); 
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
