\name{computeSensitivities}
\alias{computeSensitivities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{computeSensitivities
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeSensitivities(mathModel, pDerivsInfo, CNOlist, indexes, times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mathModel}{
%%     ~~Describe \code{mathModel} here~~
}
  \item{pDerivsInfo}{
%%     ~~Describe \code{pDerivsInfo} here~~
}
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
}
  \item{indexes}{
%%     ~~Describe \code{indexes} here~~
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
function(mathModel,pDerivsInfo,CNOlist,indexes,times)
{
  pderivs=pDerivsInfo$listPDeriv;
  numParameters=length(mathModel$index_k)+length(mathModel$index_n)+
      length(mathModel$index_k);
  numStates=length(mathModel$expressions);
  numExperiments=length(CNOlist$valueStimuli[,1]);
  simData=simulateModel(CNOlist,indexes,mathModel,times);
  numInputs=length(mathModel$inputs);
  numPderivs=numParameters*numStates;
  res=list();
  
  for(i in 1:length(mathModel$index_inh))assign(mathModel$parNames[mathModel$index_inh[i]],0);
  for(i in 1:numParameters)assign(mathModel$parNames[i],mathModel$paramValues[i]);  
  for(i in 1:numExperiments)
  { 
    expRes=vector('list',numPderivs); 
    for(j in 1:numInputs)
    {
      assign(mathModel$inputs[j],CNOlist$valueStimuli[i,j]); 
    }
    
    for(j in indexes$CNOindexInhibitors)
    {
      print(j) 
    }
    index[i]=1;
    
    for(j in 1:length(times))
    {
      for(k in 1:numStates)
      {
        assign(mathModel$outputs[k],simData[[i]][j,k]); 
      }
      count=0;
      for(k in 1:numStates)
      {
        for(r in 1:numParameters)
        {
          count=count+1;
          if(pderivs[[count]]!=0)
          {
            expRes[[count]]=c(expRes[[count]],eval(pderivs[[count]]))
          }
          else{expRes[[count]]=0;}
        }
      }
    }
    res[[i]]=expRes
  }
  return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
