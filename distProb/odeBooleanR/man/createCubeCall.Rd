\name{createCubeCall}
\alias{createCubeCall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createCubeCall(inputNames, output, truthTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputNames}{
%%     ~~Describe \code{inputNames} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
}
  \item{truthTable}{
%%     ~~Describe \code{truthTable} here~~
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
function(inputNames,output,truthTable)
{
    numInputs=length(inputNames)
    numStates=2^numInputs;
    firstFlag=TRUE;
    strODE="(";
    if(length(truthTable)>0)
    {
      for(i in 1:numStates)
      {
        if(truthTable[i])
        {
          if(firstFlag){firstFlag=FALSE;}
          else strODE=paste(strODE,"+",sep="");
          tempBin=dec2bin(i-1,numInputs)
          for(j in 1:numInputs)
          {  
            ij=tempBin[j];
            if(!as.logical(ij))
            {
              tempStr=genNormHillCubes(inputNames[j],output);
              strODE=paste(strODE,"(1-",tempStr,")",sep="");
            }
            else if(as.logical(ij))
            {
              tempStr=genNormHillCubes(inputNames[j],output);
              strODE=paste(strODE,tempStr,sep="");
            }
            if(j<numInputs) strODE=paste(strODE,"*",sep="");
          }
        }      
      }
    }
    else{strODE='0'}; 
    parTau=paste(output,"_Tau",sep="");
    parInh=paste(output,"_Inh",sep=""); 
    strODE=paste("(",strODE,"-",output,")*",parTau,")","*(1-",parInh,")",sep="");
    return(strODE);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
