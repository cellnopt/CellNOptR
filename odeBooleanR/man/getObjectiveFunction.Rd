\name{getObjectiveFunction}
\alias{getObjectiveFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getObjectiveFunction
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getObjectiveFunction(paramInfo1, y01, expData1, times1, obsIndex1, fullpars1, indexPars1, cFileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paramInfo1}{
%%     ~~Describe \code{paramInfo1} here~~
}
  \item{y01}{
%%     ~~Describe \code{y01} here~~
}
  \item{expData1}{
%%     ~~Describe \code{expData1} here~~
}
  \item{times1}{
%%     ~~Describe \code{times1} here~~
}
  \item{obsIndex1}{
%%     ~~Describe \code{obsIndex1} here~~
}
  \item{fullpars1}{
%%     ~~Describe \code{fullpars1} here~~
}
  \item{indexPars1}{
%%     ~~Describe \code{indexPars1} here~~
}
  \item{cFileName}{
%%     ~~Describe \code{cFileName} here~~
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
function(paramInfo1,y01,expData1,times1,obsIndex1,fullpars1,indexPars1,cFileName)
{
  if(cFileName=="mymod")
  {
    objectiveFunction<-function(x, y0=y01,paramInfo=paramInfo1,expData=expData1,times=times1,obsIndex=obsIndex1,fullpars=fullpars1,indexPars=indexPars1)
    {    
       fullpars[indexPars]=x;
       simData<-lsode(y0, times, func = "derivs", parms=fullpars, jacfunc=NULL, dllname="mymod",initfunc ="initmod",nout = 0);
       res=(expData[,obsIndex]-simData[,obsIndex])^2;
      
       return(res);
    }             
  }
  else if(cFileName=="myODE")
  {   
    objectiveFunction<-function(x, y0=y01,paramInfo=paramInfo1,expData=expData1,times=times1,obsIndex=obsIndex1,fullpars=fullpars1,indexPars=indexPars1)
    {    
       fullpars[indexPars]=x;
       simData<-cvodes(y0,times,"myODE","rhs",fndata =fullpars,verbose =FALSE);
       res=sum((expData[,obsIndex]-simData[,obsIndex])^2); 
       if(is.nan(res))res=10000000000;
       print(res);
       return(res);
    }
  }
  else 
  {  
    objectiveFunction<-function(x)
    {
      print("NO VALID OBJECTIVE FUNCTION WAS PROVIDED")
      return(NA);
    }
    print("A valid File Name was not provided");
  }
  
 return(objectiveFunction);   
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
