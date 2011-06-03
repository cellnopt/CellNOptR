\name{getPartialDerivs}
\alias{getPartialDerivs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getPartialDerivs
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPartialDerivs(mathModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mathModel}{
%%     ~~Describe \code{mathModel} here~~
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
function(mathModel)
{
  numParameters=length(mathModel$index_k)+length(mathModel$index_n)+
      length(mathModel$index_k);
  numStates=length(mathModel$expressions);
  listPDeriv=list();
  count=0;
  indexNonZeroPDeriv=c();
  tag=c();
  for(i in 1:numStates)
  {
    for(j in 1:numParameters)
    {
      count=count+1;
      tag=c(tag,paste('d',mathModel$outputs[i],'/d',mathModel$parNames[j],sep=""));
      listPDeriv[[count]]=D(mathModel$expressions[i],mathModel$parNames[j]);
      if(listPDeriv[[count]]!=0)indexNonZeroPDeriv=c(indexNonZeroPDeriv,count);
    }
  }
  res=list(listPDeriv=listPDeriv,tag=tag,indexNonZeroPDeriv=indexNonZeroPDeriv);
  return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
