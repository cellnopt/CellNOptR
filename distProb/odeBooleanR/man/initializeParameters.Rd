\name{initializeParameters}
\alias{initializeParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
initializeParameters(paramInfo, typeInit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paramInfo}{
%%     ~~Describe \code{paramInfo} here~~
}
  \item{typeInit}{
%%     ~~Describe \code{typeInit} here~~
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
function(paramInfo,typeInit)
{
  if(typeInit=="randomNK")
  {
    indexNK=c(paramInfo$index_k,paramInfo$index_n)
    paramInfo$paramValues[indexNK]=runif(length(indexNK));
  }
  else if(typeInit=="randomBoundedNK")
  {
    indexK=paramInfo$index_k;
    indexN=paramInfo$index_n;
    paramInfo$paramValues[indexK]=runif(length(indexK),0.5,5);
    paramInfo$paramValues[indexN]=runif(length(indexN),0.5,7);
    paramInfo$UB[indexK]=5;
    paramInfo$LB[indexK]=0.5;
    paramInfo$UB[indexN]=7;
    paramInfo$LB[indexN]=0.5;
  }
  else
  {
    print("Check you have selected any valid Option");
  }
  return(paramInfo);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
